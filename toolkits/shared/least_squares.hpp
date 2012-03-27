#ifndef GRAPHLAB_LEAST_SQUARES
#define GRAPHLAB_LEAST_SQUARES

inline void parse_edge(const edge_data& edge, const vertex_data & pdata, mat & Q, vec & vals, int i){
       
  /*for (int j=0; j<D; j++)
      Q(j,i) = pdata.pvec(j);
  }*/
  set_col(Q, i, pdata.pvec);
  vals(i) = edge.weight;  
}

float predict(const vertex_data& v1, const vertex_data& v2, float rating, float & prediction){
   //predict missing value based on dot product of movie and user iterms
   prediction = dot(v1.pvec, v2.pvec);
   //truncate prediction to allowed values
   prediction = std::min((double)prediction, mi.maxval);
   prediction = std::max((double)prediction, mi.minval);
   //return the squared error
   float sq_err = powf(prediction - rating, 2);
   return sq_err;
}


void compute_least_squares(mat & Q, vec & vals, vec & result, bool isuser, bool toprint, int numedges){
  
    //COMPUTE LEAST SQUARES (ALTERNATING LEAST SQUARES)
    //compute weighted regularization (see section 3.2 of Zhou paper)
   double reg = regularization;

   if (regnormal)
	   reg*= Q.cols();

   //std::cout<<(Q*transpose(Q)+mi.eDT*reg + mi.eDT*reg)<<std::endl;
   //std::cout<<(Q*vals)<<std::endl;
   bool ret = ls_solve(Q*transpose(Q)+ mi.eDT*reg, Q*vals, result);
   if (!ret)
     logstream(LOG_FATAL)<<"NUmerical error when computing least square step. Try to increas reglarization using --regularization=XX command line argument" << std::endl;   

  if (toprint){
    print_vec("ALS", result, false);
    std::cout << " edges: " << numedges << std::endl;
  }
}


  /***
 * UPDATE FUNCTION
 */
struct als_lapack:
  public iupdate_functor<graph_type,als_lapack>{
  void operator()(icontext_type &context){
#ifndef USE_GRAPHLAB_ENGINE
  }
  void operator+=(const als_lapack& other){ 
  }
  void finalize(iglobal_context_type& context){
  }
};

  void update_function_backslash(dummy_context& context){
#endif
  /* GET current vertex data */
  vertex_data& vdata = context.vertex_data();
 
  int id = context.vertex_id();
  bool toprint = debug && info.toprint(id);
  bool isuser = info.is_row_node(id);
  /* print statistics */
  if (toprint){
    printf("entering %s node  %u \n", (!isuser ? "movie":"user"), id);   
    debug_print_vec((isuser ? "V " : "U") , vdata.pvec, data_size);
  }

  vdata.value = 0;

  edge_list outs = context.out_edges();
  edge_list ins = context.in_edges();
  int numedges = ins.size() + outs.size();
  if (numedges == 0)
    return;

  mat Q(data_size,numedges); //linear relation matrix
  vec vals(numedges); //vector of ratings

  //USER NODES    
  if (isuser){

    for(uint i=0; i< outs.size(); i++) {
      const edge_type & edget = outs[i];
      const edge_data & edge = context.edge_data(edget);
      const vertex_data  & pdata = context.const_vertex_data(edget.target()); 
        //go over each rating of a movie and put the movie vector into the matrix Q
        //and vector vals
        parse_edge(edge, pdata, Q, vals, i); 
        if (toprint && (i==0 || i == (uint)numedges-1)){
          std::cout<<"set col: "<<i<<std::endl; 
          print_vec("", get_col(Q,i), false);
        }

      /*if (!ac.round_robin){
        gl_types::update_task task(context.target(oedgeid), user_movie_nodes_update_function);
          scheduler.add_task(task, 1);
      }*/

    }
  }

  else {


    //MOVIE NODES
    for (uint i=0; i< ins.size(); i++) {
      const edge_type & edget = ins[i];
      const edge_data & edge = context.edge_data(edget);
      const vertex_data & pdata = context.const_vertex_data(edget.source()); 
        //go over each rating by user
        parse_edge(edge, pdata, Q, vals, i); 
        if (toprint && (i==0 || i == numedges-1)){
          std::cout<<"set col: "<<i<<std::endl;
          print_vec("", get_col(Q,i), false);
        }

       float prediction;
       double trmse = predict(vdata, pdata, edge.weight, prediction); 
       if (debug && (i== (uint)info.get_start_node(false) || i == (uint)info.get_end_node(false)-1))
         cout<<"RMSE sq_err: " << trmse << " prediction: " << prediction << endl; 
  
       if (toprint)
          cout<<"trmse: " << trmse << endl;
      //aggregate RMSE
       vdata.value += trmse; 

       /*if (!ac.round_robin && trmse > ac.threshold && ps.iiter < ac.iter){
        gl_types::update_task task(context.source(iedgeid), user_movie_nodes_update_function);
          scheduler.add_task(task, 1);
      }*/
    }

  }

  vec result;
  compute_least_squares(Q, vals, result, isuser, toprint, numedges);
  vdata.pvec = result;
}
 
core<graph_type, als_lapack> * glcore_backslash = NULL;

void init_math(graph_type * _pgraph, core<graph_type, als_lapack> * _glcore, bipartite_graph_descriptor & _info, 
                   bool update_function = false){
      pgraph = _pgraph;
      glcore_backslash = _glcore;
      info = _info;
      mi.reset_offsets();
      mi.update_function = update_function;
    }


DistSlicedMat& DistSlicedMat::operator=(DistMat & mat){
			assert(mi.dist_sliced_mat_backslash);
			assert(start < end);
			assert(start_offset < end_offset);
			INITIALIZE_TRACER(als_lapack_trace, "Update function als lapack");
			BEGIN_TRACEPOINT(als_lapack_trace);
			if (mi.update_function){
				for (vertex_id_type i = start; i < (vertex_id_type)end; i++)
					glcore_backslash->schedule(i, als_lapack()); 
				runtime += glcore->start();
			}
			else {
#ifdef USE_GRAPHLAB_ENGINE
				glcore->aggregate_now("als_lapack"); 
#else
#pragma omp parallel for
				for (int t=start; t< end; t++){
					dummy_context con(t);
					update_function_backslash(con);
				}
#endif
			}
			END_TRACEPOINT(als_lapack_trace);
			//debug_print(name);
			mi.reset_offsets();
			return *this;
		}


#endif //GRAPHLAB_LEAST_SQUARES
