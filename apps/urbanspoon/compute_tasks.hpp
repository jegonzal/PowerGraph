////////////////// Initialization Task ////////////////////////////////
void init_vertex(engine_type::context_type& context,
                   graph_type::vertex_type& vertex) {
  vertex_data& vdata = vertex.data(); 
  typedef boost::unordered_map<size_t, double>::value_type kv_type;
  // initialize feature latent cache
  vdata.xfactor.setZero();
  foreach(kv_type& kv, vdata.features) {
    vdata.xfactor += kv.second * feature_latent[kv.first];
  }
  // initialize feature cache
  vdata.fcache = 0;
  typedef boost::unordered_map<size_t, double>::value_type kv_type;
  foreach(kv_type& kv, vdata.features) {
    vdata.fcache += feature_weights[kv.first] * kv.second;
  }
}

int init_edge_hold_out (graph_type::edge_type& edge) {
  // sample 20% of movies as validate set
  if (TEST_PERCENT > 0)  {
    bool flip = random::fast_bernoulli(TEST_PERCENT);
    if (flip) {
      edge.data().role = edge_data::VALIDATE;
    }
  }
  return edge.data().role == edge_data::TRAIN;
}


void init_edge(engine_type::context_type& context,
               graph_type::edge_type& edge) {
  // initialize feature cache
  edge.data().fcache = 0;
  typedef boost::unordered_map<size_t, double>::value_type kv_type;
  foreach(kv_type& kv, edge.data().features) {
    edge.data().fcache += feature_weights[kv.first] * kv.second;
  }
}

// helper function for set merge
std::set<size_t>& set_merge(std::set<size_t>& v1, const std::set<size_t>& v2) {
  v1.insert(v2.begin(), v2.end());
  return v1;
}

// Random intialize latent vector for each feature, and synchronize to all nodes. 
void init_global_vars(graphlab::distributed_control& dc) {
  w0 = 0;

  // synchronize on feature keys
  std::set<size_t> keys;
  typedef boost::unordered_map<size_t, double>::value_type kv_type;
  foreach(kv_type& kv, feature_weights) {
    keys.insert(kv.first);
  }
  dc.all_reduce2(keys, set_merge);

  // fill up feature keys
  feature_keys.clear();
  feature_keys.assign(keys.begin(), keys.end());

  // initialize feature weights and feature latent vecs.
  foreach(size_t k, feature_keys) {
    feature_latent[k].resize(vertex_data::NLATENT);
    feature_latent[k].setZero();
    feature_weights[k] = 0;
  }

  for (size_t i = 0; i < feature_keys.size(); ++i) {
    reverse_feature_keys[feature_keys[i]] = i;
  }
}


////////////////// ALS Compute Task ////////////////////////////////
/**
 * \brief The gather type used to construct XtX and Xty needed for the ALS
 * update
 *
 * To compute the ALS update we need to compute the sum of
 * \code
 *  sum: XtX = nbr.factor.transpose() * nbr.factor
 *  sum: Xy  = nbr.factor * edge.obs
 * \endcode
 * For each of the neighbors of a vertex.
 *
 * To do this in the Gather-Apply-Scatter model the gather function
 * computes and returns a pair consisting of XtX and Xy which are then
 * added. The gather type represents that tuple and provides the
 * necessary gather_type::operator+= operation.
 *
 */
class als_gather_type {
  public:
    /**
     * \brief Stores the current sum of nbr.factor.transpose() *
     * nbr.factor
     */
    mat_type XtX;

    /**
     * \brief Stores the current sum of nbr.factor * edge.obs
     */
    vec_type Xy;

    /** \brief basic default constructor */
    als_gather_type() { }

    /**
     * \brief This constructor computes XtX and Xy and stores the result
     * in XtX and Xy
     */
    als_gather_type(const vec_type& X, const double y) :
      XtX(X.size(), X.size()), Xy(X.size()) {
      XtX.triangularView<Eigen::Upper>() = X * X.transpose();
      Xy = X * y;
    } // end of constructor for gather type

    /** \brief Save the values to a binary archive */
    void save(graphlab::oarchive& arc) const { arc << XtX << Xy; }

    /** \brief Read the values from a binary archive */
    void load(graphlab::iarchive& arc) { arc >> XtX >> Xy; }

    /**
     * \brief Computes XtX += other.XtX and Xy += other.Xy updating this
     * tuples value
     */
    als_gather_type& operator+=(const als_gather_type& other) {
      if(other.Xy.size() == 0) {
        ASSERT_EQ(other.XtX.rows(), 0);
        ASSERT_EQ(other.XtX.cols(), 0);
      } else {
        if(Xy.size() == 0) {
          ASSERT_EQ(XtX.rows(), 0);
          ASSERT_EQ(XtX.cols(), 0);
          XtX = other.XtX; Xy = other.Xy;
        } else {
          XtX.triangularView<Eigen::Upper>() += other.XtX;
          Xy += other.Xy;
        }
      }
      return *this;
    } // end of operator+=
}; // end of gather type

#define ALS_MAP_REDUCE 0
/** The gather function computes XtX and Xy */
als_gather_type als_map(const graph_type::vertex_type& center,
                    graph_type::edge_type& edge,
                    const graph_type::vertex_type& other) {
  if (edge.data().role == edge_data::TRAIN) {
    vec_type other_total =  (other.data().factor + other.data().xfactor);
    double residual = (edge.data().obs - edge.data().pred);
    residual += other_total.dot(center.data().factor);
    return als_gather_type(other_total, residual);
  } else {
    return als_gather_type();
  }
} // end of gather function


als_gather_type als_map_edge (const graph_type::edge_type& e, size_t fid)  {
  vertex_data& udata = e.source().data();
  vertex_data& vdata = e.target().data();
  ASSERT_TRUE(fid != 0);
  if (e.data().role == edge_data::TRAIN && udata.features.count(fid)) {
    double residual = (e.data().obs - e.data().pred);
    double scalar = udata.features.at(fid);
    residual +=  scalar * (vdata.factor + vdata.xfactor).dot(feature_latent[fid]);
    return als_gather_type(scalar * (vdata.factor + vdata.xfactor), residual);
  } else if (e.data().role == edge_data::TRAIN && vdata.features.count(fid)) {
    double residual = (e.data().obs - e.data().pred);
    double scalar = vdata.features.at(fid);
    residual +=  scalar * (udata.factor + udata.xfactor).dot(feature_latent[fid]);
    return als_gather_type(scalar * (udata.factor + udata.xfactor), residual);
  } else  {
    return als_gather_type();
  }
}

als_gather_type als_sum(als_gather_type& v1, const als_gather_type& v2) {
  v1 += v2;
  return v1;
}

void als_update_function(engine_type::context_type& context,
                         graph_type::vertex_type& vertex,
                         const engine_type::message_type& unused) {
  // Get and reset the vertex data
  vertex_data& vdata = vertex.data();
  als_gather_type sum = context.map_reduce<als_gather_type>(ALS_MAP_REDUCE, ALL_EDGES);

  // Determine the number of neighbors.  Each vertex has only in or
  // out edges depending on which side of the graph it is located
  if(sum.Xy.size() == 0) { vdata.residual = 0; ++vdata.nupdates; return; }
  mat_type XtX = sum.XtX;
  vec_type Xy = sum.Xy;

  // Solve the least squares problem using eigen ----------------------------
    // Add regularization
  double regularization = LAMBDA;
  for(int i = 0; i < XtX.rows(); ++i)
    XtX(i,i) += regularization;
  const vec_type old_factor = vdata.factor;
  vdata.factor = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
  vdata.residual = (vdata.factor - old_factor).cwiseAbs().sum() / XtX.rows();
  ++vdata.nupdates;
}


/**************************************************************************/
/*                                                                        */
/*                     REGRESSION ON FEATURE WEIGHTS                      */
/*                                                                        */
/**************************************************************************/
als_gather_type regression_edge_map(graph_type::edge_type& edge) {
  if (edge.data().role == edge_data::TRAIN) {
    vec_type left_total = (edge.source().data().factor + edge.source().data().xfactor);
    vec_type right_total =  (edge.target().data().factor + edge.target().data().xfactor);
    double residual = (edge.data().obs - edge.data().pred);
    residual  += edge.source().data().fcache;
    residual  += edge.target().data().fcache;
    residual  += edge.data().fcache;

    vec_type x; x.resize(feature_keys.size()); x.setZero();

    typedef boost::unordered_map<size_t, double>::value_type kv_type;
    foreach(kv_type& kv, edge.source().data().features) {
      x[reverse_feature_keys[kv.first]] = kv.second;
    }
    foreach(kv_type& kv, edge.target().data().features) {
      x[reverse_feature_keys[kv.first]] = kv.second;
    }
    foreach(kv_type& kv, edge.data().features) {
      x[reverse_feature_keys[kv.first]] = kv.second;
    }

//     for (size_t i = 0; i < feature_keys.size(); ++i) {
//       size_t key = feature_keys[i];
//       if (edge.source().data().features.count(key)) {
//         x[i] = edge.source().data().features[key];
//       }
//       if (edge.target().data().features.count(key)) {
//         x[i] = edge.target().data().features[key];
//       }
//     }
    return als_gather_type(x, residual);
  } else {
    return als_gather_type();
  }
}


/**************************************************************************/
/*                                                                        */
/*                           COMPUTE BIAS TERMS                           */
/*                                                                        */
/**************************************************************************/
double compute_bias(graph_type::edge_type& edge) {
  if (edge.data().role == edge_data::TRAIN) {
    return edge.data().obs - edge.data().pred + w0;
  } else {
    return 0;
  }
}

#define BIAS_MAP_REDUCE 2
typedef std::pair<double, size_t> double_size_pair_t;
double_size_pair_t bias_map(const graph_type::vertex_type& center,
                            graph_type::edge_type& edge,
                            const graph_type::vertex_type& other) {
  if (edge.data().role == edge_data::TRAIN) {
    double res = edge.data().obs - edge.data().pred + center.data().bias;
    return double_size_pair_t (res, 1);
  } else {
    return double_size_pair_t (0, 0);
  }
}

template<typename T1, typename T2>
std::pair<T1, T2> pair_sum(std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
  v1.first += v2.first;
  v1.second += v2.second;
  return v1;
}

void compute_vertex_bias(engine_type::context_type& context,
                          graph_type::vertex_type& vertex) {
  double_size_pair_t sum = context.map_reduce<double_size_pair_t>(BIAS_MAP_REDUCE, ALL_EDGES);
  if (sum.second != 0) {
    vertex.data().bias = sum.first / sum.second;
  } else {
    vertex.data().bias = 0;
  }
}

/** Write prediction onto edges */
void compute_prediction (engine_type::context_type& context,
                         graph_type::edge_type& edge,
                         bool truncate_training = false) {
  vertex_data& udata = edge.source().data(); 
  vertex_data& vdata = edge.target().data();

  double pred = w0 + udata.bias + vdata.bias + (udata.factor + udata.xfactor)
      .dot(vdata.factor + vdata.xfactor) + vdata.fcache + udata.fcache + edge.data().fcache;

  // typedef boost::unordered_map<size_t, double>::value_type kv_type;
  // foreach(kv_type& kv, udata.features) {
  //   pred += feature_weights[kv.first] * kv.second;
  // }
  // foreach(kv_type& kv, vdata.features) {
  //   pred += feature_weights[kv.first] * kv.second;
  // }
  if (edge.data().role != edge_data::TRAIN || 
      (edge.data().role == edge_data::TRAIN && truncate_training)) {
    if (pred > MAX_VAL) {
      pred = MAX_VAL;
    } else if (pred < MIN_VAL) {
      pred = MIN_VAL;
    }
  }
  edge.data().pred = pred;
}
