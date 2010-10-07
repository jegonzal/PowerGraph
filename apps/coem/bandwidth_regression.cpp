/*
 *
 * Do regression test on pairwise bandwidth performance.
 *
 */
 

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <dirent.h>

#include <boost/random.hpp>

#include <graphlab.hpp>
 #include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
 #include <graphlab/macros_def.hpp>

using namespace std;
using namespace graphlab;
   

// Distributed control reference
int myprocid;
int numofprocs;

 

/**
  *  MAIN FUNCTION 
  */
int main(int argc,  char ** argv) {
  
  // Initialize distributed control
  graphlab::distributed_control dc(&argc, &argv);
  dc.init_message_processing(4);
  dc.barrier();
  myprocid = dc.procid();
  numofprocs = dc.numprocs();
  
  long wait_between_secs= 30;
  long duration_mins = 10;
  
  graphlab::command_line_options   clopts("Bandwidth regression test.");
  clopts.attach_option("wait_between_secs", &wait_between_secs, wait_between_secs,
  				"Wait between tests in seconds.");
  clopts.attach_option("duration", &duration_mins, duration_mins,
  				"Test duration in minutes.");
  				
 // Create a graphlab 
  if(!clopts.parse(argc, argv)) {
     std::cout << "Error in parsing input." << std::endl;
     return EXIT_FAILURE;
  }		
  
  timer t;
  t.start();
  
  char fname[255];
  sprintf(fname, "stats_bandwidth_test_%d.csv", numofprocs);
  FILE * f = fopen(fname, "w");
  
  
  while(t.current_time() < duration_mins*60) {
    dc.barrier();
    double t1 = t.current_time();
    for (size_t i = 0;i < dc.numprocs(); ++i) {
      if (dc.procid() == i)  distributed_metrics::instance(&dc)->execute_bandwith_test();
      dc.barrier();
    }
    std::cout << myprocid << ": test took " << (t.current_time()-t1) << " secs." << std::endl;
    dc.barrier();
    
    
    if (myprocid == 0) {
        for(int fromproc=0; fromproc<numofprocs; fromproc++) {
            
            for(int toproc=0; toproc<numofprocs; toproc++) {
                fprintf(f, "%lf,",  distributed_metrics::instance(&dc)->get_value(toproc, distributed_metrics::instance(&dc)->bandwidth_key(fromproc, toproc))); 
            }
            fprintf(f, "\n");  
        }
    }
    std::cout << "======== WAITING FOR NEXT TEST: " << wait_between_secs <<
        ". Test will end in " << (duration_mins*60-t.current_time())/60 << " mins" << std::endl;
    sleep(wait_between_secs);
    fflush(f);
  }  
  fclose(f);
}

 
