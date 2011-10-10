#ifndef _UNITTEST_HPP
#define _UNITTEST_HPP

extern advanced_config ac;

void test_cosamp();

/**
 * UNIT TESTING
 */
void verify_result(double obj, double train_rmse, double validation_rmse){
   assert(ac.unittest > 0);
   switch(ac.unittest){
      case 1: //ALS: Final result. Obj=0.0114447, TRAIN RMSE= 0.0033 VALIDATION RMSE= 1.1005.
	 assert(pow(obj - 0.012,2) < 1e-3);
         assert(pow(train_rmse - 0.0033,2) < 1e-3);
	 assert(pow(validation_rmse - 1.1005,2) < 1e-2);
	 break;

      case 71: //Lanczos
/*eigenvalue 0 val: 17.1428
eigenvalue 1 val: 1.83558
eigenvalue 2 val: 0.951951
eigenvalue 3 val: 0.780815
eigenvalue 4 val: 0.0329737 */
         assert(pow(get_val(ps.V,4,0) - 17.1428,2) <1e-5);
         assert(pow(get_val(ps.V,3,0) - 1.83558,2) <1e-5);
         assert(pow(get_val(ps.V,2,0) - 0.951951,2) <1e-5);
         assert(pow(get_val(ps.V,1,0) - 0.780815,2) <1e-5);
         assert(pow(get_val(ps.V,0,0) - 0.0329737,2)<1e-5);

/*   -0.1235    0.3194    0.7068   -0.3074   -0.5372
   -0.4525    0.5548   -0.4707    0.3454   -0.3830
    0.8323    0.1010   -0.0972    0.3109   -0.4370
   -0.2942   -0.6665    0.1976    0.5373   -0.3761
   -0.0279   -0.3684   -0.4800   -0.6332   -0.4819
*/       assert(ps.U.rows() == 5 && ps.U.cols() == 5);
         ASSERT_LE(norm(ps.U - init_mat("-0.1235    0.3194    0.7068   -0.3074   -0.5372 -0.4525    0.5548   -0.4707    0.3454   -0.3830 0.8323    0.1010   -0.0972    0.3109   -0.4370 -0.2942   -0.6665    0.1976    0.5373   -0.3761 -0.0279   -0.3684   -0.4800   -0.6332   -0.4819", 5, 5)) , 2);       
         break;

      case 91: //WEIGHTED_ALS: -Iter100... UV. objective=0.0207271, RMSE=0.0043/0.6344. Time to finish=0.00hr.
         // On UBUNTU 11.04 with ITPP we get: Final result. Obj=0.0154572, TRAIN RMSE= 0.0032 VALIDATION RMSE= 0.6833.
         // Ubuntu 11.04 with ITPP (Sagar Soni) we get: Final result. Obj=0.0197717, TRAIN RMSE= 0.0043 VALIDATION RMSE= 0.6678.
         //assert(pow(obj -  0.0207271,2)<1e-5);
         assert(obj < 0.021);
         //assert(pow(train_rmse - 0.0043,2)<1e-5);
         assert(train_rmse < 0.0044);
         assert(validation_rmse < 0.6840);
         //assert(pow(validation_rmse - 0.6344,2)<1e-3);
         break;

   }
}


//UNIT TESTING OF VARIOUS METHODS
void unit_testing(int unittest, command_line_options& clopts){

   if (unittest == 1){ //ALTERNATING LEAST SQUARES
      ac.datafile = "als"; ps.algorithm = ALS_MATRIX; ac.algorithm = ALS_MATRIX; ac.FLOAT = true; ac.als_lambda = 0.001;
      clopts.set_scheduler_type("round_robin(max_iterations=100,block_size=1)");
   }
   else if (unittest == 51){ //SVD++
    ac.datafile = "movielens"; ps.algorithm = SVD_PLUS_PLUS; ac.algorithm = SVD_PLUS_PLUS; clopts.set_ncpus(1); ac.debug = true; ac.maxval = 5; ac.minval = 1;
      clopts.set_scheduler_type("round_robin(max_iterations=10,block_size=1)");
  }
   else if (unittest == 71){ //Lanczos
     ac.datafile = "lanczos2"; ps.algorithm = LANCZOS;  ac.algorithm = LANCZOS; ac.matrixmarket= true; clopts.set_ncpus(1); ac.debug = true; clopts.set_scheduler_type("sweep");
   }
   else if (unittest == 91){ //WEIGHTED ALTERNATING LEAST SQUARES
      ac.datafile = "wals"; ac.algorithm = WEIGHTED_ALS; ps.algorithm = WEIGHTED_ALS; ac.FLOAT = true; ac.als_lambda = 0.001;
      clopts.set_scheduler_type("round_robin(max_iterations=100,block_size=1)");
   }
   else if (unittest == 101){
      test_cosamp();
      exit(0);
   } 
   else {
      logstream(LOG_ERROR) << " Unit test mode " << unittest << " is not supported yet! " << std::endl;
      exit(1);
   }
}



#endif
