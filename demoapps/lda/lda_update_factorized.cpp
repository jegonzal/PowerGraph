//////////////  This can't work



#include "corpus.hpp"
#include "lda.hpp"



#include <graphlab/macros_def.hpp>

extern size_t ntopics;
extern size_t nwords;
extern double;
extern double;

std::vector< glshared<size_t> > n_t;



 
class lda_update : public gl::iupdate_functor::factorized {
  size_t iters_remaining;
  std::vector<count_type> old_global_nt;
  std::vector<count_type> local_n_t;
  std::vector<count_type> local_n_dt;

public:

  lda_update(const size_t iters = 100) : iters_remaining(iters) { } 
    
  void operator+=(const lda_update& other) {
    if(!other.n_wt.empty()) {
      typedef map_type::value_type value_type;
      // Merge the other map into this map
      foreach(const value_type& pair, other.n_wt) {
        n_wt[pair.first] = pair.second;
      }
    }
  }

  
  void gather(iscope_type& scope, icallback_type& callback,
              edge_id_type eid) {
    // Make a local copy of the global topic counts
    if(local_n_t.empty()) {
      local_n_t.resize(ntopics);
      local_n_dt = scope.vertex_data()
      for(size_t t = 0; t < ntopics; ++t) 
        local_n_t[t] = global_n_t[t].get_val();
      old_global_n_t = local_n_t;
    }

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






    // only gather on a document
    ASSERT_EQ(scope.vertex_data().type == DOCUMENT);
    // Get local data structures
    const vertex_data& doc       = scope.vertex_data();
    const gl::vertex_id word_vid = scope.target(eid);
    const word_id_type word_id   = word_vid;
    const vertex_data& word      = scope.neighbor_vertex_data(word_vid);
    const edge_data& edata       = scope.edge_data(eid);

    // Initialize edge data if it is not already initialized
    if(edata.n_t.size() != ntopics) {
      edata.clear();  edata.n_t.resize(ntopics,0);
    }

    // Compute the topic probability
    std::vector<double> prob(ntopics);   
    double normalizer = 0; 
    for(size_t t = 0; t < ntopics; ++t) {
      
      const double doc_n_t = doc.n_t.size() != ntopics? 0 : doc_n_t[t];

      prob[t] = (alpha + min(doc_n_t - edata.n_t[t], 0)) * 
        (beta + min(word.n_t[t] - edata.n_t[t], 0)) / 
        (beta * nwords + min(n_t[t].get_val() - edata.n_t[t], 0) );
      normalizer += prob[t];
    }
    ASSERT_GT(normalizer, 0);
    // Normalize the probability
    for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;

    // Actually do the flipping for all occurences of the words in
    // this document
    std::vector<count_type>& new_counts = n_wt[word_id];
    new_counts.clear(); new_counts.resize(ntopics, 0);
    for(size_t i = 0; i < edata.count; ++i) {
      const size_t topic_id = graphlab::random::multinomial(prob); 
      ASSERT_LT(topic_id, ntopics);
      new_counts[topic_id]++;
    }  
  } // end of gather

  void apply(iscope_type& scope, icallback_type& callback) {
    // update the topic counts for the document
    ASSERT_EQ(scope.vertex_data().type == DOCUMENT);
    // Get local data structures
    vertex_data& doc       = scope.vertex_data();
    // Clear the document topic parameter
    doc.n_t.clear();  doc.n_t.resize(ntopics, 0);   
    // Update the document topic parameter
    typedef map_type::value_type value_type;
    foreach(const value_type& pair, other.n_wt) {
      ASSERT_EQ(pair.second.size() == ntopics);
      for(size_t t = 0; t < pair.second.size(); ++t) 
        doc.n_t[]
      n_wt[pair.first] = pair.second;
    }

   

  } // end of apply

  void scatter(iscope_type& scope, icallback_type& callback,
               edge_id_type eid) {

    // update the topic counts for the words as well as the edges 

  } // end of scatter

 
}; // end of lda_update_factorized 

