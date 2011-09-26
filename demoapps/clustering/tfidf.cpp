#include "clustering.h"
#include "../gabp/advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;



void update_weights(int * tfidf){
  for (int i=0; i<ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     for (int j=0; j< data.datapoint.nnz(); j++){
       data.datapoint.set(data.datapoint.get_nz_index(j), (data.datapoint.get_nz_data(j)/data.datapoint.nnz()) * log (ps.M/(1.0+(double)tfidf[data.datapoint.get_nz_index(j)])));
     } 
  }

 

}



void calc_tfidf(int * tfidf){
   
  memset(tfidf, 0, sizeof(double)*ps.N);
  for (int i=0; i<ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     for (int j=0; j< data.datapoint.nnz(); j++){
       tfidf[data.datapoint.get_nz_index(j)]++; 
     } 
  }

}


void tfidf_weighting(){


   logstream(LOG_INFO)<< "Computing tfidf weighting on data" << std::endl;

   int * tfidf = new int[ps.N];
   calc_tfidf(tfidf);

   update_weights(tfidf);

   delete[] tfidf;

}

