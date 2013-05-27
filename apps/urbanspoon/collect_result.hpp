/////////////////////////// Collect Top Prediction //////////////////////
#define COLLECT_TASK 2
/// Collect, Prediction, Recommendation
typedef map_join<graphlab::vertex_id_type, double> map_id_rating_t;
typedef std::pair<map_id_rating_t, map_id_rating_t> map_join_pair;
map_join_pair collect_map (const graph_type::vertex_type& center,
                           graph_type::edge_type& edge,
                           const graph_type::vertex_type& other) {
  map_join_pair ret;
  if (edge.data().role == edge_data::TRAIN) {
    ret.first.data[other.id()] = edge.data().obs; // save the old rating
  } else {
    // use prediction
    double pred = center.data().factor.dot(other.data().factor);
    ret.second.data[other.id()] = pred; // save the prediction
  }
  return ret;
}

map_join_pair collect_sum (map_join_pair& v1, const map_join_pair& v2) {
  v1.first += v2.first;
  v1.second += v2.second;
  return v1;
}

std::string rank_list_to_string(const std::vector<std::pair<double, graphlab::vertex_id_type> >& ls) {
  std::stringstream sstream;
    for(size_t i = 0;i < ls.size(); ++i) {
      graphlab::vertex_id_type gid = ls[i].second;
      sstream << "\t" << id2movieid(gid)
              << ": " << mlist[gid]
              << " = " << ls[i].first << "\n";
    }
  return sstream.str();
}

void collect_function (engine_type::context_type& context,
                       graph_type::vertex_type& vertex) {
  if (is_user(vertex)) {
    map_join_pair sum = context.map_reduce<map_join_pair>(COLLECT_TASK, ALL_EDGES);
    vertex.data().top_rated = sum.first.get_top_k(10);
    vertex.data().top_pred = sum.second.get_top_k(5);
  }
}

class recommendation_writer {
 public:
   std::string save_vertex(const graph_type::vertex_type& v) {
    std::stringstream sstream;
    if (is_user(v)) {
      const std::vector<std::pair<double, graphlab::vertex_id_type> >& top_rated 
          = v.data().top_rated;
      const std::vector<std::pair<double, graphlab::vertex_id_type> >& top_pred
          = v.data().top_pred;
      if (top_rated.size() < 10 || top_pred.size() == 0) {
        return "";
      }
      // save top rated
      sstream << v.id() << " "; 
      sstream << pair2str(top_rated[0]); 
      for (size_t i = 1; i < top_rated.size(); ++i) {
        sstream << "," <<  (pair2str(top_rated[i]));
      }
      // save top pred
      sstream << " ";
      sstream << pair2str(top_pred[0]);
      for (size_t i = 1; i < top_pred.size(); ++i) {
        sstream << "," << (pair2str(top_pred[i]));
      }
      sstream << "\n";
      return sstream.str();
    } else {
      return "";
    }
   }
   std::string save_edge(const graph_type::edge_type& e) {
     return "";
   }

 private:
   std::string pair2str(const std::pair<double, graphlab::vertex_id_type>& pair) {
     int printingid = -pair.second - SAFE_NEG_OFFSET;
     return boost::lexical_cast<std::string>(printingid) + ":" + boost::lexical_cast<std::string>(pair.first);
   }
};
