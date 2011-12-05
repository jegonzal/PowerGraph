#ifndef _UNITTEST_HPP
#define _UNITTEST_HPP

#include "clustering.h"
#include "../gabp/advanced_config.h"

extern advanced_config ac;

void test_distance();
/**
 * UNIT TESTING
 */
void verify_result(double obj, double train_rmse, double validation_rmse){
   assert(ac.unittest > 0);
   switch(ac.unittest){
      case 1: //ALS: Final result. Obj=0.0114447, TRAIN RMSE= 0.0033 VALIDATION RMSE= 1.1005.
	 break;

      case 3:
         assert(fabs(ps.cost - 0.522652) < 1e-5);
         break;


      case 101:
         //three clusters, 0->19, 20-> 39, 40 -> 59,
	 assert(get_val(ps.output_assignements,0,0) == 0);
         assert(get_val(ps.output_assignements,19,0)== 0);
      	 assert(get_val(ps.output_assignements,20,0) == 2);
         assert(get_val(ps.output_assignements,39,0) == 2);
       	 assert(get_val(ps.output_assignements,40,0) == 1);
         assert(get_val(ps.output_assignements,59,0) == 1);
	 break;

      case 102:
         //three clusters, 0->19, 20-> 39, 40 -> 59,
	 assert(get_val(ps.output_assignements,0,0) > 0.9);
         assert(get_val(ps.output_assignements,19,0) > 0.9);
      	 assert(get_val(ps.output_assignements,20,1) > 0.9);
         assert(get_val(ps.output_assignements,39,1) > 0.9);
       	 assert(get_val(ps.output_assignements,40,2) > 0.9);
         assert(get_val(ps.output_assignements,59,2) > 0.9);
	 break;
         
   }
}


//UNIT TESTING OF VARIOUS METHODS
void unit_testing(int unittest, graphlab::command_line_options& clopts){
   if (unittest == 1){
      test_math();
      exit(0);
   }
   else if (unittest == 2){
      test_distance();
      exit(0);
   }
   /*if (unittest == 1){ //ALTERNATING LEAST SQUARES
      //ac.datafile = "als"; ps.algorithm = ALS_MATRIX; ac.algorithm = ALS_MATRIX; ac.FLOAT = true; ac.als_lambda = 0.001;
      //clopts.set_scheduler_type("round_robin(max_iterations=100,block_size=1)");
   }*/
   else if (unittest == 3){
      ac.datafile = "cluster5x5";
      ac.algorithm = K_MEANS_FUZZY;
      ac.iter = 10;
      ac.init_mode = 0;
      ac.debug=true;
      clopts.set_ncpus(1);
      ac.K = 3;
   
   }
   else if (unittest == 4){
     test_fmath();
     exit(0);
   }
   else if (unittest ==70){
     ac.datafile = "lanczos2"; ac.algorithm = SVD_EXPERIMENTAL; ac.K = 2; ac.init_mode = 0; ac.matrixmarket = true;
     ac.debug =true; 
   }
   else if (unittest ==71){
     ac.datafile = "lanczos2"; ac.algorithm = SVD_EXPERIMENTAL; ac.K=2; ac.init_mode = 0; ac.matrixmarket = true; ac.svd_step=1; ac.debug=true;
   }
   else if (unittest == 72){
     ac.datafile = "lanczos2"; ac.algorithm = SVD_EXPERIMENTAL; ac.K=2; ac.init_mode = 0; ac.matrixmarket = true; ac.svd_step=2; ac.debug=true;
   }
   else if (unittest == 50){
//r netflix 5 2 0 --pmfformat=true --float=false --ncpus=8 --knn_sample_percent=0.8
     ac.datafile = "netflix"; ac.algorithm = ITEM_KNN; ac.K = 2; ac.init_mode =0; ac.supportgraphlabcf = true; ac.isfloat=false; ac.ncpus=8; ac.knn_sample_percent=0.8;
     clopts.set_ncpus(8);
   }
   else if (unittest == 51){
     ac.datafile = "netflix"; ac.algorithm = USER_KNN; ac.K = 2; ac.init_mode =0; ac.supportgraphlabcf = true; ac.isfloat=false; ac.ncpus=8; ac.knn_sample_percent=0.02;
     clopts.set_ncpus(8);
   }
   else if (unittest == 71){
//./glcluster  lanczos2 7 2 0 --matrixmarket=true --svd_compile_eigenvectors=true  --reduce_mem_consumption=true --debug=true --svd_compile_eigenvectors_block_size=1 
     ac.datafile = "lanczos2"; ac.algorithm = SVD_EXPERIMENTAL; ac.K = 2; ac.init_mode = 0; ac.matrixmarket = true;
     ac.svd_compile_eigenvectors_block_size = 1; ac.svd_compile_eigenvectors = true; ac.reduce_mem_consumption = true; ac.debug = true;
   }
  else if (unittest == 101){
//./glcluster clive_test 0 3 1 --pmfformat=true --float=true --distance_metric=8
      ac.datafile = "clive_test";
      ac.algorithm = K_MEANS;
      ac.K = 3;
      ac.init_mode = 1;
      ac.supportgraphlabcf = true;
      ac.FLOAT = true;
      ac.distance_measure = 8;
   }
   else if (unittest == 102){
 //clive_test 2 3 0 --pmfformat=true --float=true --fuzzy_exponent=1.1 --debug=true --max_iter=10 --ncpus=1
      ac.datafile = "clive_test";
      ac.algorithm = K_MEANS_FUZZY;
      ac.K = 3;
      ac.init_mode = 0;
      ac.supportgraphlabcf = true;
      ac.FLOAT = true;
      ac.fuzzy_exponent = 1.3;

   }
   else {
      logstream(LOG_ERROR) << " Unit test mode " << unittest << " is not supported yet! " << std::endl;
      exit(1);
   }
 }



#endif
