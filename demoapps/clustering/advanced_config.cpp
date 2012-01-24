#include "clustering.h"
#include "distance.h"
#include "../gabp/advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;
extern const char * inittypenames[];
extern const char * distance_measure_name[];

void advanced_config::init_command_line_options(graphlab::command_line_options & clopts){

  clopts.attach_option("datafile", &datafile, datafile, "Input matrix/tensor");
  clopts.add_positional("datafile");

  clopts.attach_option("algorithm", &algorithm, algorithm, "algorithm: 0=K-means, 1=Kmeans++, 2=Fuzzy K-mean, 3=LDA 4=K-Shell decomposition, 5=Item-KNN, 6=User-KNN, 7 = SVD (experimental)");
  clopts.add_positional("algorithm");
  
  clopts.attach_option("K", &K, K, "number of clusters");
  clopts.add_positional("K");  

  clopts.attach_option("init_mode", &init_mode, init_mode, "Initialization: 0 = random, 1= round_robin, 2 = kmeans++, 3=random cluster");
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
  clopts.attach_option("svd_step", &svd_step, svd_step, "svd_step. 0 = full computation, 1 = iterations only, 2 = eigendecomposition only"); 
  clopts.attach_option("svd_orthogonolize", &svd_orthogonolize, svd_orthogonolize, "perform orthogonalization in SVD");
  clopts.attach_option("fuzzy_exponent", &fuzzy_exponent, fuzzy_exponent, "Fuzzy K-means exponent (between 1.01 -> 2)");
  clopts.attach_option("fuzzy_scatter", &fuzzy_scatter, fuzzy_scatter, "Fuzzy K-means random initialization variance (between 0.01 -> 2)");
  clopts.attach_option("init_clusters_from_file", &init_clusters_from_file, init_clusters_from_file, "K-MEANS: Init cluster heads from file");
  clopts.attach_option("N", &N, N, "feature width size (number of matrix columns");
}

void problem_setup::verify_setup(){

   if (ac.datafile.size() == 0){
      logstream(LOG_ERROR) << "Input file is not specified. Aborting run. " << std::endl;
      exit(1);
   }

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
  if (ac.distance_measure != EUCLIDEAN)
    logstream(LOG_INFO) << "Setting distance metric to : " << distance_measure_name[ac.distance_measure] << std::endl;

  if (ac.init_mode == INIT_KMEANS_PLUS_PLUS && ac.algorithm != K_MEANS_PLUS_PLUS){
     logstream(LOG_WARNING) << "When running with --init_mode=2 (Kmeans++), algorithm should be Kmeans++ (--algorithm=1)" << std::endl;
     ac.algorithm = K_MEANS_PLUS_PLUS;
  }

  output_assignements_comment = "%%GraphLab clustering library output file\n";
  output_clusters_comment = output_assignements_comment;

  switch(algorithm){
    case USER_KNN:
    case ITEM_KNN:
      output_assignements_comment += std::string("%%This file includes the ") + boost::lexical_cast<std::string>(ac.K) + std::string(" closest points to each user/item where the points are numbered between 1 to max_point\n%%First column = data point number, Second Column is the index 1..K, Third column is the closest data point\n%%Note that -1 assignemnt means the point had no non-zero feature and thus not assigned to any other point\n");
      output_clusters_comment += std::string("%%This file includes the distances (using ") + std::string(distance_measure_name[ac.distance_measure]) +
		                 std::string(") to each point in the assignemnets file\n%%First column is data point number, second column is the index 1..K \
	, third column is the distance\n%%Note that -1 means that this point was not assigned to any other close point since it has all zero features\n");
      break;

   case K_MEANS_PLUS_PLUS:
   case K_MEANS:
      output_assignements_comment += std::string("%%This file contains cluster assignment for each of the data points\n");
      output_clusters_comment += std::string("%%This file contains clusters center for each of the ") + boost::lexical_cast<std::string>(ac.K) + std::string("clusters\n");
      break;


   case K_MEANS_FUZZY:
      output_clusters_comment += std::string( "%%This file contains ") + boost::lexical_cast<std::string>(ac.K) + std::string(" cluster locations for fuzzy k-means\n");
      output_assignements_comment += std::string("%%This file contains weighted assignemtns of each data point to each fizzy cluster\n");
      break;


   case SVD_EXPERIMENTAL:
     output_comment3 = output_assignements_comment;
     output_comment4 = output_assignements_comment;
     V_comment += std::string("%%This file contains the matrix V which is received by [U,D,V']=svd(A)\n");
     U_comment += std::string("%%This file contains the matrix U which is received by [U,D,V']=svd(A)\n");
     output_comment3 += std::string("%%This file contains the eigenvectors of AA'\n");
     output_comment4 += std::string("%%This file contains the eigenvectors of A'A\n");
      break;

   default:
      break;
   };



  switch(algorithm){
     case K_MEANS_PLUS_PLUS:
     case K_MEANS:
     case ITEM_KNN:
     case USER_KNN:
       output_assignements_integer = true; // is the output integer or float?
       break;

     default:
       output_assignements_integer = false;

  } 

  if (algorithm == SVD_EXPERIMENTAL){
    if (ac.svd_step == SVD_FULL_PASS){
      ac.reduce_mem_consumption = false;
    }
    else if (ac.svd_step == SVD_ITERATIONS){
      ac.reduce_mem_consumption = true;
      ac.svd_finalize = false;
      ac.svd_compile_eigenvectors = false;
    }
    else if (ac.svd_step == SVD_EIGENDECOMPOTISION){
      ac.reduce_mem_consumption = true;
      ac.svd_finalize = true;
      ac.svd_compile_eigenvectors = true;
    }
  }


  if (ac.init_clusters_from_file && ac.algorithm == K_MEANS){
    ps.init_type = INIT_FROM_FILE;
    assert(ac.N > 0);
    ps.N = ac.N;
  }
}
