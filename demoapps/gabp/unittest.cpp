#include "linear.h"
#include "advanced_config.h"

extern advanced_config config;
extern graphlab::glshared<int> MATRIX_WIDTH_KEY;
extern problem_setup ps;

void unittest(graphlab::command_line_options &clopts){
  if (config.unittest == 1){
    logstream(LOG_WARNING)<< "Going to run GaBP unit testing using matrix of size 3x2" << std::endl;
    //const char * args[] = {"gabp", "0", "mat3x2", "1e-10", "--unittest=1","--syncinterval=10000"};
    //clopts.parse(6, (char**)args);
    config.algorithm = GaBP;
    config.datafile = "mat3x2";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    //config.oldformat=true;
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 2){
    logstream(LOG_WARNING)<< "Going to run GaBP unit testing using matrix of size 3x3" << std::endl;
    //const char * args[] = {"gabp", "0", "mat3x3", "1e-10", "--unittest=2", "--syncinterval=10000"};
    //clopts.parse(6, (char**)args);
    config.algorithm = GaBP;
    config.datafile = "mat3x3";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    //config.oldformat = true;
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 3){
    logstream(LOG_WARNING)<< "Going to run Jacobi unit testing using matrix of size 3x3" << std::endl;
    //const char * args[] = {"gabp", "1", "mat3x3", "1e-10", "--unittest=3", "--syncinterval=10000"};
    //clopts.parse(6, (char**)args);
    config.algorithm = JACOBI;
    config.datafile = "mat3x3";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    //config.oldformat = true;
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 4){
    logstream(LOG_WARNING)<< "Going to run CG unit testing using matrix of size 3x3" << std::endl;
    //const char * args[] = {"gabp", "2", "mat3x3", "1e-10", "--unittest=4", "--syncinterval=10000"};
    //clopts.parse(6, (char**)args);
    config.iter = 10;
    config.algorithm = CONJUGATE_GRADIENT;
    config.datafile = "mat3x3";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    //config.oldformat = true;
    clopts.set_scheduler_type("fifo");
  }
 else if (config.unittest == 5){
    logstream(LOG_WARNING)<< "Going to run CG unit testing using matrix of size 3x2" << std::endl;
    //const char * args[] = {"gabp", "2", "mat3x2", "1e-10", "--unittest=5", "--syncinterval=10000"};
    config.algorithm = CONJUGATE_GRADIENT;
    config.datafile = "mat3x2";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    // clopts.parse(6, (char**)args);
    config.iter = 10;
    //config.oldformat = true;
    clopts.set_scheduler_type("fifo");
  }
 else if (config.unittest == 6){
    logstream(LOG_WARNING)<< "Going to run GaBP inverse unit testing using matrix of size 3x3" << std::endl;
    //const char * args[] = {"gabp", "3", "mat3x3", "1e-10", "--unittest=6", "--syncinterval=10000"};
    //clopts.parse(6, (char**)args);
    config.algorithm = GaBP_INV;
    config.datafile = "mat3x3";
    config.threshold = 1e-10;
    config.syncinterval = 10000;
    config.iter = 10;
    //config.oldformat = true;
    MATRIX_WIDTH_KEY.set(3);
    clopts.set_scheduler_type("round_robin(max_iterations=20,block_size=1)");
 }
 else if (config.unittest == 7){
     logstream(LOG_WARNING)<< "Going to run Shotgun lasso with arcene dataset" << std::endl;
    config.algorithm = SHOTGUN_LASSO;
    config.datafile = "arcene";
    config.threshold = 1e-5;
    config.iter = 10;
    config.display_cost = true;
    clopts.set_ncpus(1);
    config.shotgun_lambda = 1;
 }
 else if (config.unittest == 8){
     logstream(LOG_WARNING)<< "Going to run Shotgun lasso with arcene dataset" << std::endl;
    config.algorithm = SHOTGUN_LOGREG;
    config.datafile = "arcene";
    config.threshold = 1e-5;
    config.iter = 10;
    config.display_cost = true;
    clopts.set_ncpus(1);
    config.shotgun_lambda = 1;
 }


}

void verify_unittest_result(double diff){
   if (config.unittest == 1){
      assert(diff <= 1.7e-4);
   }
   else if (config.unittest == 2){
      assert(diff <= 1e-15);
   }
   else if (config.unittest == 3){
      assert(diff <= 1e-30);
   }
   else if (config.unittest == 4 || config.unittest == 5)
      assert(diff < 1e-14);
   else if (config.unittest == 7){
      assert(pow(ps.last_cost - 97.3638, 2)< 1e-6);
   }
   else if (config.unittest == 8){
      assert(pow(ps.last_cost - -65.1957, 2) <1e-6);
   } 
}
