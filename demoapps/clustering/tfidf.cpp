#include "clustering.h"
#include "../gabp/advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;



void update_weights(int * tfidf){
  for (int i=0; i<ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     FOR_ITERATOR(j, data.datapoint){
       //data.datapoint.set(data.datapoint.get_nz_index(j), (data.datapoint.get_nz_data(j)/data.datapoint.nnz()) * log (ps.M/(1.0+(double)tfidf[data.datapoint.get_nz_index(j)])));
       set_div(data.datapoint, j, 1.0/nnz(data.datapoint)* log (ps.M/(1.0+(double)tfidf[get_nz_index(data.datapoint, j)])));
     } 
  }

 

}



void calc_tfidf(int * tfidf){
   
  memset(tfidf, 0, sizeof(double)*ps.N);
  for (int i=0; i<ps.M; i++){
     vertex_data & data = ps.g<graph_type>()->vertex_data(i);
     FOR_ITERATOR(j, data.datapoint){
       tfidf[get_nz_index(data.datapoint, j)]++; 
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

