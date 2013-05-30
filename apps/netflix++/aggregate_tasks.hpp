//////////////////////////// RMSE Aggregator /////////////////////////////
struct error_aggregator : public graphlab::IS_POD_TYPE {
  double train;
  double test;
  size_t ntrain;
  size_t ntest;
  error_aggregator() : train(0), test(0), ntrain(0), ntest(0) {}
  error_aggregator& operator+=(const error_aggregator& other) {
    train += other.train;
    test += other.test;
    ntrain += other.ntrain;
    ntest += other.ntest;
    return *this;
  }
};

/**
 * \brief Given an edge compute the error associated with that edge
 */
error_aggregator extract_l2_error(const graph_type::edge_type & edge) {
  double err = (edge.data().obs - edge.data().pred) * (edge.data().obs - edge.data().pred);
  error_aggregator ret;
  if (edge.data().role == edge_data::TRAIN) {
    ret.train = err;
    ret.ntrain = 1;
  } else {
    ret.test = err;
    ret.ntest = 1;
  }
  return ret;
} // end of extract_l2_error

//////////////////////////// Norm Residual Aggregator /////////////////////////////
double extract_residual(const graph_type::vertex_type& vertex) {
  return vertex.data().residual;
}
double extract_l2_norm(const graph_type::vertex_type& vertex) {
  return vertex.data().factor.norm() + vertex.data().xfactor.norm();
}

//////////////////////////// Recommendation Aggregator /////////////////////////////
template<typename key_t, typename val_t>
struct map_join {
  boost::unordered_map<key_t, val_t> data; // rating than movie
  void save(graphlab::oarchive& oarc) const {
    oarc << data;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> data;
  }
  map_join& operator+=(const map_join& other) {
    std::copy(other.data.begin(), other.data.end(),
              std::inserter(data, data.end()));
    return *this;
  }

  template<typename Comparator>
  std::vector<std::pair<key_t, val_t> > get_top_k(size_t n, Comparator cmp) {
    std::vector<std::pair<key_t, val_t> > ret;
    typename boost::unordered_map<key_t, val_t>::const_iterator iter = data.begin();
    std::vector<std::pair<key_t, val_t> > all_copy;
    while (iter != data.end()) {
      all_copy.push_back(std::make_pair(iter->first, iter->second));
      ++iter;
    }
    std::sort(all_copy.rbegin(), all_copy.rend(), cmp); 
    size_t limit = all_copy.size() < n ? all_copy.size() : n;
    std::copy(all_copy.begin(), all_copy.begin() + limit,
              std::inserter(ret, ret.end()));
    return ret;
  }

  template<typename Comparator>
  std::vector<std::pair<key_t, val_t> > get_top_k_exclude(
      size_t n, std::set<key_t>& exclude, Comparator cmp) {
    std::vector<std::pair<key_t, val_t> > ret;
    typename boost::unordered_map<key_t, val_t>::const_iterator iter = data.begin();
    std::vector<std::pair<key_t, val_t> > all_copy;
    while (iter != data.end()) {
      if (exclude.count(iter->first) == 0) {
	      all_copy.push_back(std::make_pair(iter->first, iter->second));
      }
      ++iter;
    }
    std::sort(all_copy.rbegin(), all_copy.rend(), cmp);
    size_t limit = all_copy.size() < n ? all_copy.size() : n;
    std::copy(all_copy.begin(), all_copy.begin() + limit,
              std::inserter(ret, ret.end()));
    return ret;
  }
};
