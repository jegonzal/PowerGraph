// /////////////////////// Collect Explanation for Top Predictions //////////////
// #define EXPLAIN_TASK 3
// /// Collect Items to be Explained
// struct rec_items {
//   std::vector<vec_type> xs; 
//   std::vector<vertex_id_type> ids;
// 
//   size_t size() { return ids.size(); }
// 
//   rec_items() { }
//   rec_items (const vec_type& v, const vertex_id_type id)  {
//     xs.push_back(v);
//     ids.push_back(id);
//   }
// 
//   /** \brief Save the values to a binary archive */
//   void save(graphlab::oarchive& arc) const { 
//     arc << xs << ids; 
//   }
// 
//   /** \brief Read the values from a binary archive */
//   void load(graphlab::iarchive& arc) { 
//     arc >> xs >> ids; 
//   }
// };
// 
// rec_items gather_rec_items (const graph_type::vertex_type& center,
//                     graph_type::edge_type& edge,
//                     const graph_type::vertex_type& other) {
//   ASSERT_TRUE(is_user(center));
//   if (edge.data().role == edge_data::prediction) {
//     bool recommended = false;
//     // check whether the movie is recommended in topK 
//     for (size_t i = 0; i < center.top_pred.size(); ++i) {
//       if (other.id() == center.top_pred[i].second) {
//         recommended = true;
//         break;
//       }
//     }
//     if (recommended) {
//       return rec_items(other.data().factor, other.id());
//     }
//   }
//   return rec_items();
// }
// 
// rec_items combine_rec_items (rec_items& v1, const rec_items& v2) {
//   v1.xs.insert(v1.xs.end(), v2.xs.begin(), v2.xs.end());
//   v1.ids.insert(v1.ids.end(), v2.ids.begin(), v2.ids.end());
//   return v1;
// }
// 
// /// Collect explanation to the items
// typedef map_join<vertex_id_type, double> map_id_score_t;
// 
// // vector of size == #target_items, with element i stores topK explanatory 
// // items for target_items[i]
// typedef std::vector<map_id_score_t> explain_gather_t;
// 
// explain_gather_t
//   gather_explain_score(const graph_type::vertex_type& center,
//                      graph_type::edge_type& edge,
//                      const graph_type::vertex_type& other) {
//   ASSERT_TRUE(is_user(center));
//   explain_gather_t ret;
//   if (edge.data().role ==  edge_data::TRAIN) {
//     // compute explantory score for each target item
//     vec_type scores = center.data().xWu * other.data().factor;
//     for (size_t i = 0; i < scores.size(); ++i) {
//       map_id_score_t m;
//       m.data[other.id()] = scores[i];
//       ret.push_back(m);
//     }
//   }
//   return ret;
// }
// 
// explain_gather_t
//   combine_explain_score(explain_gather_t& v1, const explain_gather_t& v2) {
//     if (v2.size() == 0) {
//       return v1;
//     }
//     if (v1.size() == 0) {
//       return v2;
//     }
//     ASSERT_EQ(v1.size() == v2.size());
//     for (size_t i = 0; i < v1.size(); ++i) {
//       v1[i] += v2[i];
//     }
//     return v1;
// }

struct factor_gather_t {
  map_join<vertex_id_type, vec_type> pred;
  map_join<vertex_id_type, std::pair<vec_type, double> > train;

  void save(oarchive& arc) const {
    arc << pred << train;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(iarchive& arc) {
    arc >> pred >> train; 
  }
};


#define FACTOR_GATHER_TASK 3
factor_gather_t factor_gather(const graph_type::vertex_type& center,
                                    graph_type::edge_type& edge,
                                    const graph_type::vertex_type& other) {
  ASSERT_TRUE(is_user(center));
  factor_gather_t ret;
  if (edge.data().role == edge_data::TRAIN) {
    ret.train.data[other.id()] = std::make_pair<vec_type, double> (
        other.data().factor, edge.data().obs);
  } else {
    ret.pred.data[other.id()] = other.data().factor;
  }
  return ret;
}

factor_gather_t factor_combine(factor_gather_t& v1, const factor_gather_t& v2) {
  v1.pred += v2.pred;
  v1.train += v2.train;
  return v1;
}

void exp_collect_function (engine_type::context_type& context,
                         graph_type::vertex_type& vertex) {
  if (!is_user(vertex)) { return; }
  vertex_data& vdata = vertex.data();
  vdata.top_explain.clear();

  if (vdata.top_pred.size() == 0) { return; }

  factor_gather_t factors = context.map_reduce<factor_gather_t>(FACTOR_GATHER_TASK, ALL_EDGES);
  mat_type Wu = vdata.Wu.selfadjointView<Eigen::Upper>();
  Wu = Wu.inverse();

  // for each recommended movie
  for (size_t i = 0; i < vdata.top_pred.size(); ++i) {
    map_join<vertex_id_type, double> scores;
    const vec_type& xj = factors.pred.data[vdata.top_pred[i].second];
    // compute contribution of each rated movie
    boost::unordered_map<vertex_id_type, std::pair<vec_type, double> >::const_iterator it = factors.train.data.begin();
    while (it != factors.train.data.end()) {
      vertex_id_type id = it->first;
      const vec_type& xi = it->second.first;
      double rating = it->second.second;
      ASSERT_EQ(xi.size(), xj.size());
      double score = (xi.transpose() * Wu).dot(xj) * rating;
      scores.data[id] = score;
      ++it;
    }
    // record the 3 most important movies
    vdata.top_explain.push_back(scores.get_top_k(3));
  }
}

std::string explain_list_to_string(
    const std::vector< std::pair<double, graphlab::vertex_id_type> >& rec_ls,
    const std::vector< std::vector<std::pair<double, graphlab::vertex_id_type> > >& exp_ls) {
  ASSERT_EQ(rec_ls.size(), exp_ls.size());
  std::stringstream sstream;
    for(size_t i = 0; i < rec_ls.size(); ++i) {
      graphlab::vertex_id_type gid = rec_ls[i].second;
      sstream << id2movieid(gid) << " " << mlist[gid] << ":\n";
      for (size_t j = 0; j < exp_ls[i].size(); ++j) {
        graphlab::vertex_id_type gid2 = exp_ls[i][j].second;
        double score = exp_ls[i][j].first;
        sstream << "\t" << id2movieid(gid2) << ":" << mlist[gid2] << "=" << score <<"\n";
      }
    }
  return sstream.str();
}

class explain_writer {
 public:
   std::string save_vertex(const graph_type::vertex_type& v) {
     std::stringstream sstream;
     if (is_user(v)) {
       if (v.data().top_rated.size() < 10 || v.data().top_pred.size() == 0) {
         return "";
       }
       const std::vector<std::vector<std::pair<double, graphlab::vertex_id_type> > >& top_explain = v.data().top_explain;
       const std::vector<std::pair<double, graphlab::vertex_id_type> >& top_pred= v.data().top_pred;
       std::stringstream sstream;
       for(size_t i = 0; i < top_pred.size(); ++i) {
         graphlab::vertex_id_type gid = top_pred[i].second;
         for (size_t j = 0; j < top_explain[i].size(); ++j) {
           graphlab::vertex_id_type gid2 = top_explain[i][j].second;
           double score = top_explain[i][j].first;
           sstream << id2movieid(gid) << "\t" << id2movieid(gid2) << "\t" << score <<"\n";
         }
       }
       return sstream.str();
     } else {
       return "";
     }
   }

   std::string save_edge(const graph_type::edge_type& e) {
     return "";
   }
};
