#ifndef _TENSOR_HPP
#define _TENSOR_HPP

extern advanced_config ac;
extern problem_setup ps;



template<typename graph_type, typename edge_data>
const edge_data * find_edge(int time, const graph_type*g, edge_id_t edge){
     return &g->edge_data(edge);
}

template<>
const edge_data_mcmc * find_edge<>(int time, const graph_type_mult_edge*g, edge_id_t edge){
//find the right edge which matches the current time
     const multiple_edges * medges= &g->edge_data(edge);
      const edge_data_mcmc *data = NULL;
      for (int j=0; j< (int)medges->medges.size(); j++){
        data = &medges->medges[j];
        if (data->time == time)
          return data;
      }
     return NULL;
}


 
#include <graphlab/macros_def.hpp>
std::vector<edge_id_t> * edges;  //edge list for pointing from time nodes into their connected users/movies (for tensor)
double muT = 1; //mean of time nodes

// Calculate time nodes (for tensor)
template<typename graph_type, typename vertex_data, typename edge_data>
void calc_T(int i){

  assert(ps.tensor);
 
  assert(i >=0 && i < ps.K);
  if (ac.zero && edges[i].size() == 0){
     if (i == ps.K-1)
        last_iter<graph_type, vertex_data, edge_data>();
    return;
  }

  assert(edges[i].size() > 0);

  if (ac.debug && (i==0 || i == ps.K-1))
    printf("node %d with Q size: %d %d\n", i, (int)edges[i].size(), (int)ac.D);


  const graph_type * g = ps.g<graph_type>(TRAINING);

  int k=0;
  int batchSize = 100;
  mat Q(batchSize, ac.D); 
  vec vals(batchSize);
  mat QQ = zeros(ac.D,ac.D);
  vec RQ = zeros(ac.D);
  int cnt =0;


  timer t; t.start();
  foreach (edge_id_t edge, edges[i]){
    if (k % batchSize == 0){
      Q = zeros(batchSize, ac.D);
      vals = zeros(batchSize);
      cnt = 1;
    }

   const edge_data *data = find_edge<graph_type, edge_data>(i,g, edge);
   assert(data != NULL);
   assert(data->time == i); 

    assert((int)g->target(edge)>= ps.M);
    assert((int)g->source(edge)< ps.M);
    const vertex_data * v1 = &g->vertex_data(g->target(edge));
    const vertex_data * v2 = &g->vertex_data(g->source(edge));
    vec ret = elem_mult(v1->pvec, v2->pvec);
     
    for (int s=0; s<ac.D; s++)
      Q(k%batchSize,s)=ret(s);
    if (ac.debug && (i==0 || i == ps.K-1) && (k == 0 || k == (int)edges[i].size() - 1))
      std::cout<<" clmn "<<k<< " vec: " << ret<<std::endl;

    vals[k%batchSize] = data->weight;
    k++;

    if ((cnt  == batchSize) || (cnt < batchSize && k == (int)edges[i].size()-1)){
      QQ += transpose(Q)*Q;
      RQ += transpose(Q)*vals;
      assert(QQ.rows() == ac.D && QQ.cols() == ac.D);
    }
    cnt++;
  }
  ps.counter[BPTF_TIME_EDGES] += t.current_time();


  if (ac.debug && (i == 0 ||i== ps.K-1 )){
    std::cout<<"QQ:"<<QQ<<std::endl;
    std::cout<<"RQ:"<<RQ<<std::endl;
  }

  assert(RQ.size() == ac.D);
  assert((unsigned int)k == edges[i].size());
  vec sol;    
  vec out;

  t.start();
  if (i == 0){
    vec t1 = ps.times[1].pvec; 
    //calc least squares estimation of time nodes
    if (!ps.BPTF){
      QQ = QQ+ 2*ps.eDT;
      bool ret = ls_solve(QQ, ps.pT*(t1 + ps.vones*muT)+ RQ, out);
      assert(ret);
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = 2*A_T+alpha*QQ;
      mat iAk_;
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T*(t1 + mu_T)+alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, ac.D);
    }
  }
  else if (i == ps.K-1){
    vec tk_2 = ps.times[ps.K-2].pvec; 
    //calc least squares estimation of time nodes
    if (!ps.BPTF){
      QQ = QQ + ps.eDT;
      bool ret = ls_solve(QQ, ps.pT*tk_2 + RQ, out);
      assert(ret); 
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = A_T+alpha*QQ;
      mat iAk_; 
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T*tk_2+alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, ac.D);
    }
  }
  else {
    vec tsum; 
    vec t1 = ps.times[i-1].pvec; 
    vec t2 = ps.times[i+1].pvec;
    tsum = t1+t2;
    //calc least squares estimation of time nodes
    if (!ps.BPTF){
      QQ = QQ + 2*ps.eDT;
      bool ret = ls_solve(QQ, ps.pT*tsum + RQ, out);
      assert(ret);
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = 2*A_T+alpha*QQ;
      mat iAk_;
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T* tsum +alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, ac.D);
    }
  }

  ps.times[i].pvec = out;
  if (ac.debug && (i == 0|| i == ps.K-1)) 
    std::cout<<ps.times[i].pvec<<std::endl;
  assert(QQ.rows() == ac.D && QQ.cols() == ac.D);
  ps.counter[BPTF_LEAST_SQUARES] += t.current_time();
  //}

  if (i == ps.K-1){
    last_iter<graph_type, vertex_data, edge_data>();
  }
          


}
// update function for time nodes
// this function is called only in tensor mode
void time_node_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler) {

  assert(ps.tensor);

  int id = scope.vertex() - ps.M-ps.N;
  assert(id >=0 && id < ps.K);	
  if (ac.debug && (id == 0 || id == ps.K-1)){
    printf("entering time node  %d \n", id);   
  } 
  if (ps.K > 1)
    calc_T<graph_type, vertex_data, edge_data>(id); 
  else last_iter<graph_type, vertex_data, edge_data>();
}
// update function for time nodes
// this function is called only in tensor mode
void time_node_update_function(gl_types_mcmc::iscope &scope, gl_types_mcmc::icallback &scheduler) {

  assert(ps.tensor);

  int id = scope.vertex() - ps.M-ps.N;
  assert(id >=0 && id < ps.K);	
  if (ac.debug && (id == 0 || id == ps.K-1)){
    printf("entering time node  %d \n", id);   
  } 
  if (ps.K > 1)
    calc_T<graph_type_mcmc, vertex_data, edge_data_mcmc>(id); 
  else last_iter<graph_type_mcmc, vertex_data, edge_data_mcmc>();
}

// update function for time nodes
// this function is called only in tensor mode
void time_node_update_function(gl_types_mult_edge::iscope &scope, gl_types_mult_edge::icallback &scheduler) {
  assert(ps.tensor);

  int id = scope.vertex() - ps.M-ps.N;
  assert(id >=0 && id < ps.K);	
  if (ac.debug && (id == 0 || id == ps.K-1)){
    printf("entering time node  %d \n", id);   
  } 
  if (ps.K > 1)
    calc_T<graph_type_mult_edge, vertex_data, edge_data_mcmc>(id); 
  else last_iter<graph_type_mult_edge, vertex_data, multiple_edges>();
}

// update function for time nodes
// this function is called only in tensor mode
void time_node_update_function(gl_types_svdpp::iscope &scope, gl_types_svdpp::icallback &scheduler) {
  assert(ps.tensor);
}





#include <graphlab/macros_undef.hpp>
#endif //_TENSOR_HPP
