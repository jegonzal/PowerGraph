#ifndef __SVD_HPP
#define __SVD_HPP


float LRATE1 = 0.003;               // Learning rate parameter
float LRATE2 = 0.003;               // Learning rate parameter for biases
float LRATE3 = 0.001;               // Learning rate parameter for weights
float LAMBDA1 = 0.015;               // Regularization parameter used to minimize over-fitting
float LAMDA2 = 0.005;               // Biases regularization
float LAMDA3 = 0.015;               // Regularization parameter for weights

void calc_user_moviebag(vertex_data & vdata, edge_list & outs){

#ifdef GL_SVD_PP
    vdata.weight = zeros(D);
    foreach(graphlab::edge_id_t oedgeid, outs) {
       //edge_data & edge = scope.edge_data(oedgeid);
       const vertex_data  & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
       vdata.weight += pdata.weight;
    }
#endif //GL_SVD_PP
} 



float predict_svd_rating(vertex_data & user, vertex_data & movie){

  double sum = itpp::sum(itpp.dot(movie.pvec, user.pvec + itpp.dot(user.weight, user.ei)));

  sum += user.bias + movie.bias;

  //todo: sum += cashe_rate_residuals

  //todo: truncate values!
}


/***
 * UPDATE FUNCTION
 */
void svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  
  /* print statistics */
  if (debug&& (scope.vertex() == 0 || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1) || ((int)scope.vertex() == 93712))){
    printf("SVDPP: entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

  user.rmse = 0;

  int numedges = vdata.num_edges;
  if (numedges == 0){
    return; //if this user/movie have no ratings do nothing
  }


  edge_list outs = scope.out_edge_ids();
  edge_list ins = scope.in_edge_ids();
  timer t;

  int i=0;

  t.start(); 
  //USER NODES    
  if ((int)scope.vertex() < M){

    calc_user_moviebag(vdata, outs);

    foreach(graphlab::edge_id_t oedgeid, outs) {

#ifndef GL_NO_MULT_EDGES
      multiple_edges &medges =scope.edge_data(oedgeid);
#else
      edge_data & edge = scope.edge_data(oedgeid);
#endif
      const vertex_data  & movie = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
#ifndef GL_NO_MULT_EDGES	   
      for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#endif
        //go over each rating of a movie and put the movie vector into the matrix Q
        //and vector vals
        float p = predict_svd_rating(vdata, pdata);    
        float err = edge.weight -  p;
        user.rmse += err*err;

        float cf_bias = user.bias;
        float mf_bias = movie.bias;

        user.bias += (float)(LRATE2 * (err-LAMBDA2 * cf_bias));
        movie.bias += (float)(LRATE2 * (err-LAMBDA2 * mf_bias));

	//cache off old feature values
        vec cf = user.pvec;
        vec mf = movie.pvec;
        vec wf = movie.weight;

        float ei = 1.0 / sqrtf(outs.size() + 1.0); //regularization

        user.pvec += (LRATE1 * (err*mf - LAMBDA1*cf));
        movie.pvec += (LRATE1 * (err*(itpp::dot(cf+ei,user.weight)) - LAMBDA1*mf));
        movie.weight += LRATE3*(itpp::dot(err*user.ei, mf)-LAMBDA3*wf);
        user.weight += movie.weight - wf;
	i++;
#ifndef GL_NO_MULT_EDGES
      }   
#endif
           
            
    }
  }

  assert(i == numedges);
  counter[EDGE_TRAVERSAL] += t.current_time();

  calc_users_moviebag();


  vec result;
      
  if (!BPTF){
    //COMPUTE LEAST SQUARES (ALTERNATING LEAST SQUARES)
    t.start();
    double regularization = LAMBDA;
    //compute weighted regularization (see section 3.2 of Zhou paper)
    if (!regnormal)
	regularization*= Q.cols();

    bool ret = itpp::ls_solve(Q*itpp::transpose(Q)+eDT*regularization, Q*vals, result);
    assert(ret);
    counter[ALS_LEAST_SQUARES] += t.current_time();
  }
  else {
    //COMPUTE LEAST SQUARES (BPTF)
    //according to equation A.6 or A.7 in Xiong paper.
    assert(Q.rows() == D);
    t.start();
    mat iAi_;
    bool ret =inv(((int)scope.vertex() < M? A_U : A_V) + alpha *  Q*itpp::transpose(Q), iAi_);
    assert(ret);
    t.start();
    vec mui_ =  iAi_*((((int)scope.vertex() < M)? (A_U*mu_U) : (A_V*mu_V)) + alpha * Q * vals);
    counter[BPTF_LEAST_SQUARES2]+= t.current_time();
       
    t.start();
    result = mvnrndex(mui_, iAi_, D); 
    assert(result.size() == D);
    counter[BPTF_MVN_RNDEX] += t.current_time();
  }

  if (debug && (((int)scope.vertex()  == 0) || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1))){
    std::cout <<(BPTF?"BPTF":"ALS")<<" Q: " << Q << std::endl<<" result: " << result << " edges: " << numedges << std::endl;
  }
      
  //store new result
  vdata.pvec =  result;

  //calc post round tasks
  if (!tensor && (int)scope.vertex() == M+N-1)
    last_iter();

}

void last_iter(){
  printf("Entering last iter with %d\n", iiter);

  double res,res2;
  double rmse = calc_rmse_q(res);
  //rmse=0;
  printf("%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", gt.current_time(), BPTF?"BPTF":"ALS", iiter,calc_obj(res),  rmse, calc_rmse(&validation_graph, true, res2));
  iiter++;
  if (iiter == BURN_IN && BPTF){
    printf("Finished burn-in period. starting to aggregate samples\n");
  }
         
  if (BPTF){
    timer t;
    t.start();
    if (iiter > BURN_IN)
    	sample_alpha(res);
    sample_U();
    sample_V();
    if (tensor) 
      sample_T();
    counter[BPTF_SAMPLE_STEP] += t.current_time();
    if (infile == "kddcup" || infile == "kddcup2")
	export_kdd_format(&test_graph, false);
   }
}




#endif //__SVD_HPP
