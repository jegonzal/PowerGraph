struct query_result {
  vertex_data vdata;
  edge_data edata;
  void save(oarchive& oarc) const {
    oarc << vdata << edata;
  }
  void load(iarchive& iarc) {
    iarc >> vdata >> edata;
  }
};

bool compare_by_obs(const std::pair<graphlab::vertex_id_type, query_result>& v1, 
                    const std::pair<graphlab::vertex_id_type, query_result>& v2) {
  return v1.second.edata.obs < v2.second.edata.obs;
}

bool compare_by_pred(const std::pair<graphlab::vertex_id_type, query_result>& v1,
                     const std::pair<graphlab::vertex_id_type, query_result>& v2) {
  return v1.second.edata.pred < v2.second.edata.pred;
}

bool select_training(const edge_data& edata) {
  return edata.role == edge_data::TRAIN;
}

bool select_testing(const edge_data& edata) {
  return !select_training(edata);
}

map_join<graphlab::vertex_id_type, query_result> gather_neighbors(
    distributed_control& dc,
    graph_type& graph,
    graphlab::vertex_id_type& uid,
    boost::function<bool(const edge_data&)> edge_select_fn) {

  // gather all neighbors
  map_join<graphlab::vertex_id_type, query_result> all_neighbors;
  if (graph.contains_vertex(uid)) {
    graph_type::vertex_type vtx(graph.vertex(uid));
    graph_type::local_vertex_type lvtx(vtx);
    foreach(graph_type::local_edge_type edge, lvtx.out_edges()) {
      if (edge_select_fn(edge.data())) {
        graph_type::local_vertex_type target = edge.target();
        graph_type::vertex_id_type gid = target.global_id();
        query_result q;
        q.vdata = target.data();
        q.edata = edge.data();
        all_neighbors.data[gid] = q; 
      }
    }
  }
  dc.all_reduce(all_neighbors);
  return all_neighbors;
}

map_join<graphlab::vertex_id_type, query_result> query(distributed_control& dc,
           graph_type& graph,
           graphlab::vertex_id_type uid,
           size_t topk)  {
  // Gather training edges
  map_join<graphlab::vertex_id_type, query_result> training_set = gather_neighbors(dc, graph, uid, select_training);

  // Gather testing edges
  map_join<graphlab::vertex_id_type, query_result> testing_set = gather_neighbors(dc, graph, uid, select_testing);

  // output training
  std::vector<std::pair<graphlab::vertex_id_type, query_result> > topk_training = training_set.get_top_k(topk, compare_by_obs);
  for (size_t i = 0; i < topk_training.size(); ++i) {
    query_result& q = topk_training[i].second;
    size_t gid = vid2movieid(topk_training[i].first); 
        dc.cout() << gid << "\t"
        << mv_names[gid] << "\t"
        << q.edata.obs << "\t"
        << q.edata.pred << "\t"
        << 1 // is_training
        << "\n"; 
  }
  dc.cout() << "\n";

  // output testing 
  std::vector<std::pair<graphlab::vertex_id_type, query_result> > topk_testing = testing_set.get_top_k(topk, compare_by_pred);
  for (size_t i = 0; i < topk_testing.size(); ++i) {
    query_result& q = topk_testing[i].second;
    graphlab::vertex_id_type gid = vid2movieid(topk_testing[i].first); 
        dc.cout() << gid << "\t"
        << mv_names[gid] << "\t"
        << q.edata.obs  << "\t"
        << q.edata.pred << "\t"
        << 0 // not is_training
        << "\n";
  }
  dc.cout() << "\n";
  return training_set;
}

std::vector<std::pair<graphlab::vertex_id_type, double> > get_recommendation(
                        distributed_control& dc,
                        graph_type& graph, 
                        graphlab::vertex_id_type uid,
                        size_t topk,
                        std::set<graphlab::vertex_id_type>& exclude) {

  map_join<graphlab::vertex_id_type, double> all_predict;
  // now for the recommendations.
  bool is_master = graph.contains_vertex(uid) &&
      graph_type::local_vertex_type(graph.vertex(uid)).owned();

  // broadcast the user vector
  vec_type factor;
  if (is_master) {
    factor = graph.vertex(uid).data().factor;
  }
  dc.broadcast(factor, is_master);

  // now loop through all the vertices.
  for (size_t i = 0;i < graph.num_local_vertices(); ++i) {
    graph_type::local_vertex_type lvtx(graph.l_vertex(i));
    if (lvtx.owned() && (int)(lvtx.global_id()) < 0) {
      double pred = lvtx.data().factor.dot(factor);
      all_predict.data[lvtx.global_id()] = pred;
    }
  }
  dc.all_reduce(all_predict);

  std::vector<std::pair<graphlab::vertex_id_type, double> > topk_rec = 
      all_predict.get_top_k_exclude(topk, exclude, 
                                    (boost::bind(&std::pair<graphlab::vertex_id_type, double>::second, _1) < boost::bind(&std::pair<graphlab::vertex_id_type, double>::second, _2)));

  dc.cout() << "Top recommendations: \n"; 
  for (size_t i = 0; i < topk_rec.size(); ++i) {
    vertex_id_type gid = vid2movieid(topk_rec[i].first); 
    dc.cout() << gid << "\t"
              << mv_names[gid] << "\t"
              << topk_rec[i].second
              << "\n";
  }
  dc.cout() << "\n";
  return topk_rec;
}

struct explanation {
  graphlab::vertex_id_type src, dest;
  double obs, score;
  bool operator<(const explanation& e2) const {
    if (obs == e2.obs) {
      return score < e2.score;
    } else {
      return obs < e2.obs;
    }
  }
};

// Returns the explanation of vid using given training set.
std::vector<explanation> get_explanation(distributed_control& dc,
                                         graph_type& graph,
                                         graphlab::vertex_id_type vid,
                                         map_join<graphlab::vertex_id_type, query_result>& training_set) {
  als_gather_type acc;
  typedef boost::unordered_map<graphlab::vertex_id_type, query_result>::value_type kv_type;

  foreach (kv_type& kv, training_set.data) {
    acc += als_gather_type(kv.second.vdata.factor, kv.second.edata.obs);
  }

  for(int i = 0; i < acc.XtX.rows(); ++i)
    acc.XtX(i,i) += LAMBDA;
  mat_type Winv = acc.XtX.selfadjointView<Eigen::Upper>();
  Winv = Winv.inverse();

  // now for the recommendations.
  bool is_master = graph.contains_vertex(vid) &&
      graph_type::local_vertex_type(graph.vertex(vid)).owned();
  // broadcast the user vector
  vec_type factor;
  if (is_master) {
    factor = graph.vertex(vid).data().factor;
  }
  dc.broadcast(factor, is_master);

  std::vector<explanation> explain_scores;
  foreach (kv_type& kv, training_set.data) {
   double score = (factor.transpose() * Winv * kv.second.vdata.factor)[0] * kv.second.edata.obs;
   explanation exp;
   exp.src = kv.first;
   exp.dest = vid; 
   exp.obs = kv.second.edata.obs;
   exp.score = score;
   explain_scores.push_back(exp);
  }
  std::sort(explain_scores.rbegin(), explain_scores.rend());

  dc.cout() << "Top explanation to : " << mv_names[vid2movieid(vid)] << "\n"; 
  for (size_t i = 0; i < std::min(explain_scores.size(), (size_t)5); ++i) {
    graphlab::vertex_id_type gid = vid2movieid(explain_scores[i].src); 
        dc.cout() << gid << "\t"
        << mv_names[gid] << "\t"
        << explain_scores[i].obs << "\t"
        << explain_scores[i].score << "\n";
  }
  dc.cout() << "\n";

  return explain_scores;
}
