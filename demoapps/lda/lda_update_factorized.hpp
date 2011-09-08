#include "corpus.hpp"
#include "lda.hpp"



#include <graphlab/macros_def.hpp>


 
class lda_update : public gl::iupdate_functor {
  size_t iters_remaining;

  std::vector<count_type> old_global_n_t;
  std::vector<count_type> local_n_t;
  

public:

  lda_update(const size_t iters = 100) : iters_remaining(iters) { } 
    
  void operator+=(const lda_update& other) { }

  bool is_factorizable() const { return true; }

  gl::iupdate_functor::edge_set gather_edges() { 
    return gl::iupdate_functor::OUT_EDGES;
  }
  
  bool writable_gather() const { return true; }

  
  gl::iupdate_functor::edge_set scatter_edges() { 
    return gl::iupdate_functor::NO_EDGES;
  }


  void gather(iscope_type& scope, icallback_type& callback,
              edge_id_type eid) {

    if(local_n_t.empty()) {
      local_n_t.resize(ntopics);
      old_global_n_t.resize(ntopics);
      for(size_t t = 0; t < ntopics; ++t) 
        old_global_n_t[t] = local_n_t[t] = global_n_t[t].get_val();
    }
    

    // Get the data ---------------------------------------------------------
    // Get local data structures
    vertex_data& doc       = scope.vertex_data();
    // only gather on a document
    ASSERT_EQ(doc.type(), DOCUMENT);

    const gl::vertex_id word_vid = scope.target(eid);
    const word_id_type word_id   = word_vid;
    vertex_data& word      = scope.neighbor_vertex_data(word_vid);
    edge_data& edata       = scope.edge_data(eid);
    edata.init(ntopics); // ensure that the edge data is initialized
    ASSERT_LT(word_id, nwords);      
    ASSERT_EQ(word.type(), WORD);

    std::vector<double> prob(ntopics); 
    std::vector<count_type> n_dwt(ntopics);
    double normalizer = 0; 


    // Compute the probability table ----------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {
      // number of tokens from document d with topic t
      const double n_dt = 
        std::max(doc.get(t) - edata.get(t), 0);
      // number of token with word w and topic t
      const double n_wt = 
        std::max(word.get(t) - edata.get(t), 0);
      // number of tokens with topic t
      const double n_t = 
        std::max(local_n_t[t] - edata.get(t), 0);
      // compute the final probability
      prob[t] = (alpha + n_dt) * (beta + n_wt) / 
        (beta * nwords + n_t);
      normalizer += prob[t];
    }
    ASSERT_GT(normalizer, 0);
    // Normalize the probability
    for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;
      
    // Draw new topic assignments -------------------------------------------
    n_dwt.clear(); n_dwt.resize(ntopics, 0);
    for(count_type i = 0; i < edata.get_count(); ++i) {
      const size_t topic_id = graphlab::random::multinomial(prob); 
      ASSERT_LT(topic_id, ntopics);
      n_dwt[topic_id]++;
    }

    // Update all the tables ------------------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {
      doc.add(t, n_dwt[t] - edata.get(t));
      word.add(t, n_dwt[t] - edata.get(t));
      local_n_t[t] += (n_dwt[t] - edata.get(t));
      edata.set(t, n_dwt[t]);
    }

    // Reschedule self if necessary
    if(iters_remaining > 0) 
      callback.schedule(scope.vertex(), lda_update(iters_remaining-1));

  } // end of gather


  
  void apply(iscope_type& scope, icallback_type& callback) { 
    ASSERT_EQ(global_n_t.size(), ntopics);
    if(!local_n_t.empty()) {
      // update the global variables
      for(size_t t = 0; t < ntopics; ++t) {
        global_n_t[t] += (local_n_t[t] - old_global_n_t[t]);
      }    
    }
  }

  // void scatter(iscope_type& scope, icallback_type& callback,
  //              edge_id_type eid) { } // end of scatter

 
}; // end of lda_update_factorized 



#include <graphlab/macros_undef.hpp>


