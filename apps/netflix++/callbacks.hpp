#include "json_util.hpp"

class user_query_obj;

// call backs returns the list of movie names
std::pair<std::string, std::string>
movie_list_callback(std::map<std::string, std::string>& varmap) {
  const std::pair<std::string, std::string>
    pair("text/html", map2json(mv_names));
  return pair;
}

// Wrapper for rpc gather neighbors
map_join<graphlab::vertex_id_type, query_result>
  gather_training_edges_wrapper (graphlab::vertex_id_type& uid) {
    return gather_local_neighbors(uid, select_training);
}
map_join<graphlab::vertex_id_type, query_result>
  gather_testing_edges_wrapper (graphlab::vertex_id_type& uid) {
    return gather_local_neighbors(uid, select_testing);
}

std::string query_result2json (size_t vid, const query_result& q) {
  std::stringstream ss;
  size_t id = vid2movieid(vid);
  ss << "{\"id\": " << "\"" << id << "\","
     << "\"name\": " << "\"" << mv_names[id] << "\", "
     << "\"rating\": " << "\"" << q.edata.obs << "\", "
     << "\"pred\": " << "\"" << q.edata.pred << "\", "
     << "\"training\": " << (q.edata.role==edge_data::TRAIN) << "}";
  return ss.str();
}

/**
 * Return json string of a user result including
 * training/testing/recommendation/explanation 
 */
std::pair<std::string, std::string>
user_query_callback(std::map<std::string, std::string>& varmap) {
  std::pair<std::string, std::string> ret("text/html", "");
  if (varmap.count("uid") == 0) {
    ret.second = "Couldn't find uid field.";
  } else {
    // do some thing
    graphlab::vertex_id_type uid = boost::lexical_cast<graphlab::vertex_id_type>(varmap["uid"]);

    // gather training set
    map_join<graphlab::vertex_id_type, query_result> training_set = gather_training_edges_wrapper(uid);
    for (size_t i = 1; i < dc_ptr->numprocs(); ++i) {
      training_set += dc_ptr->remote_request(i, gather_training_edges_wrapper, uid);
    }
    // gather testing set
    map_join<graphlab::vertex_id_type, query_result> testing_set = gather_testing_edges_wrapper(uid);
    for (size_t i = 1; i < dc_ptr->numprocs(); ++i) {
      testing_set += dc_ptr->remote_request(i, gather_testing_edges_wrapper, uid);
    }

    // get top k
    std::vector<std::pair<graphlab::vertex_id_type, query_result> > topk_training = training_set.get_top_k(TOPK, compare_by_obs);
    std::vector<std::pair<graphlab::vertex_id_type, query_result> > topk_testing = testing_set.get_top_k(TOPK, compare_by_pred);

    // build json string
    std::vector<std::string> tmp; 
    for (size_t i = 0; i < topk_training.size(); ++i) {
      tmp.push_back(query_result2json(topk_training[i].first, topk_training[i].second));
    }
    for (size_t i = 0; i < topk_testing.size(); ++i) {
      tmp.push_back(query_result2json(topk_testing[i].first, topk_testing[i].second));
    }
    ret.second = "[" + boost::algorithm::join(tmp, ",") + "]\n";

    logstream(LOG_INFO) << "Done query uid " <<  uid << std::endl;
  }
  return ret;
}

void launch_custom_metric_server() {
  graphlab::launch_metric_server();
  graphlab::add_metric_server_callback("movie_list", movie_list_callback);
  graphlab::add_metric_server_callback("user_query", user_query_callback);
  // graphlab::add_metric_server_callback("ldaparam", lda::set_param_callback);
  // graphlab::add_metric_server_callback("lockword", lda::lock_word_callback);
  // graphlab::add_metric_server_callback("wordclouds", lda::word_cloud_callback);
  // graphlab::add_metric_server_callback("addtopic", lda::add_topic_callback);
}

// class user_query_obj {
//  public:
//    user_query_obj(distributed_control& dc): rpc(dc, this) {}
// 
//    /**
//     * RPC function synchronously called on all machines to complete the query
//     */
//    void user_query_exec(graphlab::vertex_id_type uid) {
//      rpc.barrier();
//      std::cout << "proc " << rpc.procid() << " enters barrier\n";
//      std::vector<int> values;
//      values.resize(rpc.numprocs());
//      values[rpc.procid()] = rpc.procid();
//      rpc.all_gather(values);
//      std::string ret;
//      for (size_t i = 0; i < values.size(); ++i) {
//        ret += boost::lexical_cast<std::string>(values[i]);
//      }
//      result_string = ret;
//      rpc.barrier();
//      std::cout << "proc " << rpc.procid() << " exits barrier\n";
//    }
// 
//    std::string get_result_str() {
//      return result_string;
//    }
//  private:
//    mutable dc_dist_object<user_query_obj> rpc;
//    std::string result_string;
// };
// 
// void user_query_rpc(size_t uid) {
//   user_query_obj_ptr->user_query_exec(uid);
// }
// 
