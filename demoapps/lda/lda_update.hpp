#include "corpus.hpp"
#include "lda.hpp"



#include <graphlab/macros_def.hpp>




class lda_update : public gl::iupdate_functor {  
  size_t iters_remaining;

public:

  lda_update(const size_t iters = 100) : iters_remaining(iters) { }
  

  void operator()(gl::iscope& scope, gl::icallback& callback) {
    std::cout << "Processing vertex: " << scope.vertex() << std::endl;
    ASSERT_GT(iters_remaining, 0);
    // Make a local copy of the global topic counts
    std::vector<count_type> local_n_t(ntopics);
    for(size_t t = 0; t < ntopics; ++t) 
      local_n_t[t] = global_n_t[t].get_val();
    const std::vector<count_type> old_global_n_t = local_n_t;
    

    // Get local data structures
    vertex_data& doc       = scope.vertex_data();
    // only gather on a document
    ASSERT_EQ(doc.type(), DOCUMENT);
    // Loop over the words in the document (encoded by out edges)
    const gl::edge_list out_edges = scope.out_edge_ids();
    std::vector<double> prob(ntopics); 
    std::vector<count_type> n_dwt(ntopics);
    double normalizer = 0; 
    foreach(gl::edge_id eid, out_edges) {
      // Get the data ---------------------------------------------------------
      const gl::vertex_id word_vid = scope.target(eid);
      const word_id_type word_id   = word_vid;
      vertex_data& word      = scope.neighbor_vertex_data(word_vid);
      edge_data& edata       = scope.edge_data(eid);
      edata.init(ntopics); // ensure that the edge data is initialized
      ASSERT_LT(word_id, nwords);      
      ASSERT_EQ(word.type(), WORD);

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
    } // end of for loop

    // update the global variables
    for(size_t t = 0; t < ntopics; ++t) {
      global_n_t[t] += (local_n_t[t] - old_global_n_t[t]);
    }

    // Reschedule self if necessary
    if(--iters_remaining > 0) 
      callback.schedule(scope.vertex(), *this);
    
  } // end of operator()
 
}; // end of lda_update
#include <graphlab/macros_undef.hpp>


