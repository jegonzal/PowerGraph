#include "clustering.h"

#include "../gabp/advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;
extern const char * inittypenames[];

void advanced_config::init_command_line_options(graphlab::command_line_options & clopts){

  clopts.attach_option("datafile", &datafile, datafile, "Input matrix/tensor");
  clopts.add_positional("datafile");

  clopts.attach_option("algorithm", &algorithm, algorithm, "algorithm: 0=K-means, 1=Kmeans++, 2=Fuzzy K-mean, 3=LDA");
  clopts.add_positional("algorithm");
  
  clopts.attach_option("K", &K, K, "number of clusters");
  clopts.add_positional("K");  

  clopts.attach_option("init_mode", &init_mode, init_mode, "Initialization: 0 = random, 1= round_robin, 2 = kmeans++");
  clopts.add_positional("init_mode");
  clopts.attach_option("max_iter", &iter, iter, "Maximum number if iterations");
  clopts.attach_option("debug", &debug, debug, "Display debug output. (optional)");
  clopts.attach_option("zero", &zero, zero, "support zero edges");  
  clopts.attach_option("savegraph", &savegraph, savegraph, "save graphs to file");  
  clopts.attach_option("loadgraph", &loadgraph, loadgraph, "load graphs to file");  
  clopts.attach_option("savefactors", &savefactors, savefactors, "save factors to file");  
  clopts.attach_option("loadfactors", &loadfactors, loadfactors, "load factors from file");  
  clopts.attach_option("stats", &stats, stats, "compute graph statistics");  
  clopts.attach_option("binaryoutput", &binaryoutput, binaryoutput, "export clusters and assignments into a binary file"); 
  clopts.attach_option("matrixmarket", &matrixmarket, matrixmarket, "give input in matrix market format"); 
 
  clopts.attach_option("scalerating", &scalerating, scalerating, "scale data by this factor");  
  clopts.attach_option("maxval", &maxval, maxval, "maximal allowed value in data");
  clopts.attach_option("minval", &minval, minval, "minimal allowed value in data");
  clopts.attach_option("unittest", &unittest, unittest, "unit testing. ");
  clopts.attach_option("pmfformat", &supportgraphlabcf, supportgraphlabcf, "input file of graphlab collaborative filtering library ");
  clopts.attach_option("float", &FLOAT, FLOAT, "data in float format (pmf data)"); 
  clopts.attach_option("printhighprecision", &printhighprecision, printhighprecision, "print RMSE output with high precision");
  clopts.attach_option("clusterdump", &clusterdump, clusterdump, "dump cluster locations into a text file");
  clopts.attach_option("tfidf", &tfidf, tfidf, "tfidf weighting scheme applied to document/term matrix");
  clopts.attach_option("omp_support", &omp_support, omp_support, "OMP parallelization support");
  clopts.attach_option("lda_inner_em_iter", &em_max_inner_iter, em_max_inner_iter, "LDA inner EM max iterations");
  clopts.attach_option("threshold", &threshold, threshold, "Convergence threshold");
  clopts.attach_option("distance_metric", &distance_measure, distance_measure, "distance metric");
  clopts.attach_option("knn_sample_percent", &knn_sample_percent, knn_sample_percent, "knn sample percentage (0 -> 1)");
  clopts.attach_option("reduce_mem_consumption", &reduce_mem_consumption, reduce_mem_consumption, "reduce memory consumption (potentially slower run)");
}

void problem_setup::verify_setup(){
   K = ac.K;
   algorithm = (runmodes)ac.algorithm;
   logstream(LOG_INFO) << "Setting cluster initialization mode to: " << inittypenames[ac.init_mode] << std::endl;
   ps.init_type = (initizliation_type)ac.init_mode;
   assert(K>0);

   if (ac.algorithm == K_MEANS_FUZZY){
      if (ps.init_type != INIT_RANDOM_CLUSTER){
	 logstream(LOG_WARNING) << "With Fuzzy k-means, currently supported init method is init_random cluster" << std::endl;
          ps.init_type = INIT_RANDOM_CLUSTER;
      }
   }
  if (ac.knn_sample_percent != 1.0){
     if (ps.algorithm != USER_KNN && ps.algorithm != ITEM_KNN){
        logstream(LOG_WARNING) << "Knn sample percent has no effect" << std::endl;
     }
       
  }
}
