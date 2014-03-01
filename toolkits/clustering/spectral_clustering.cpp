/*
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <graphlab.hpp>

#include <graphlab/util/fs_util.hpp>
#include <boost/filesystem.hpp>

//remove assigned options from arguments
std::string get_arg_str_without(int argc, char** argv,
    std::vector<std::string> remove_opts) {
  std::stringstream strm;
  bool skip_next = false;
  for (int i = 1; i < argc; ++i) {
    bool skip = false;
    for (size_t j = 0; j < remove_opts.size(); ++j) {
      std::string with_equal = remove_opts[j] + "=";
      if (strncmp(with_equal.c_str(), argv[i], with_equal.size()) == 0) {
        skip = true;
      } else if (strncmp(remove_opts[j].c_str(), argv[i], remove_opts[j].size())
          == 0) {
        skip = true;
        skip_next = true;
      }
    }
    if (skip == false && skip_next == false) {
      strm << argv[i] << " ";
    } else if (skip == false && skip_next == true) {
      skip_next = false;
    }
  }
  return strm.str();
}

bool call_graph_laplacian_construction(const std::string& mpi_args,
    const std::string& filename, const float sigma, const float epsilon,
    const size_t num_nearests, const std::string& args) {
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << "./graph_laplacian_for_sc ";
  strm << " --data=" << filename;
  strm << " --sigma=" << sigma;
  strm << " --similarity-thres=" << epsilon;
  strm << " --t-nearest=" << num_nearests;
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }

  return true;
}

void make_initial_vector_file(const std::string& filename, const size_t num_data){
  std::ofstream ofs((filename + ".init").c_str());
  for(size_t i=0;i<num_data;++i){
    ofs << 0.1*((i+1)%10)/10.0 << "\n";
  }
  ofs.close();
}

bool call_svd(const std::string& mpi_args, const std::string& filename,
    const std::string& svd_dir, const size_t num_clusters, const size_t rank,
    const size_t num_data, const std::string& args) {
  make_initial_vector_file(filename, num_data+1);
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << svd_dir << "svd " + filename + ".glap";
  strm << " --rows=" << num_data+1;
  strm << " --cols=" << num_data;
  strm << " --nsv=" << num_clusters;
  strm << " --nv=" << rank;
//  strm << " --tol=1e-10";
//  strm << " --max_iter=20";
  strm << " --quiet=1";
  strm << " --save_vectors=1";
  strm << " --ortho_repeats=3";
  strm << " --id=1";
  strm << " --prediction=" << filename << ".";
  strm << " --initial_vector=" << filename + ".init";
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }

  return true;
}

bool call_eigen_vector_normalization(const std::string& mpi_args,
    const std::string& filename, const std::string& graph_analytics_dir,
    const size_t num_clusters, const size_t rank, const size_t num_data,
    const std::string& args) {
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << graph_analytics_dir << "eigen_vector_normalization";
  strm << " --data=" << filename;
  strm << " --clusters=" << num_clusters;
  strm << " --rank=" << rank;
  strm << " --data-num=" << num_data;
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }

  return true;
}

bool call_kmeans(const std::string& mpi_args, const std::string& filename,
    const std::string& kmeans_dir, const size_t num_clusters,
    const std::string& args) {
  //call svd
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << kmeans_dir << "kmeans ";
  strm << " --data " << filename << ".compressed";
  strm << " --clusters " << num_clusters;
  strm << " --output-data " << filename << ".result";
  strm << " --id=1";
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }
  return true;
}


bool call_kmeans_as_preprocess(const std::string& mpi_args, const std::string& filename,
    const std::string& kmeans_dir, const size_t num_clusters,
    const std::string& args) {
  //call svd
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << kmeans_dir << "kmeans ";
  strm << " --data " << filename;
  strm << " --clusters " << num_clusters;
  strm << " --output-data " << filename << ".pre.labels";
  strm << " --output-clusters " << filename << ".pre.centers";
  strm << " --id=1";
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }
  return true;
}


//select good value of rank (TODO)
int get_lanczos_rank(const size_t num_clusters, const size_t num_data) {
  size_t rank = 1;
  if (num_data < 1000) {
    if (num_clusters + 10 <= num_data)
      rank = num_clusters + 10;
    else
      rank = num_data;
  } else if (num_data < 10000) {
    rank = num_clusters + 100;
  } else if (num_data < 100000) {
    rank = num_clusters + 150;
  } else if (num_data < 1000000) {
    rank = num_clusters + 200;
  } else {
    rank = num_clusters + 300;
  }
  return rank;
//  return num_clusters+2;
}

void read_pairs_with_prefix(std::vector<std::vector<size_t> >& ret, const std::string& prefix){
  std::string directory_name;
  std::string original_path(prefix);
  boost::filesystem::path path(prefix);
  std::string search_prefix;
  if (boost::filesystem::is_directory(path)) {
    // if this is a directory
    // force a "/" at the end of the path
    // make sure to check that the path is non-empty. (you do not
    // want to make the empty path "" the root path "/" )
    directory_name = path.native();
  }
  else {
    directory_name = path.parent_path().native();
    search_prefix = path.filename().native();
    directory_name = (directory_name.empty() ? "." : directory_name);
  }
  std::vector<std::string> files;
  graphlab::fs_util::list_files_with_prefix(directory_name, search_prefix, files);
  if (files.size() == 0) {
    logstream(LOG_WARNING) << "No files found matching " << original_path << std::endl;
  }
  for(size_t i = 0; i < files.size(); ++i) {
    std::ifstream ifs(files[i].c_str());
    if (!ifs) {
      std::cout << "can't read " << files[i] << std::endl;
      return;
    }
    while( !ifs.eof() ) {
      std::vector<size_t> pair;
      size_t id = 0;
      size_t label = 0;
      ifs >> id;
//      ifs.ignore(1);
      ifs >> label;
      if(id > 0 && label > 0){
        pair.push_back(id);
        pair.push_back(label);
        ret.push_back(pair);
      }
    }
  }
}

int recover_labels(const std::string& prefix){
  const std::string kmeans_result_prefix = prefix + ".pre.labels";
  const std::string spectral_result_prefix = prefix + ".pre.centers.result";
  const std::string outfile = prefix + ".result_1_of_1";

  std::vector<std::vector<size_t> > kmeans_result;
  read_pairs_with_prefix(kmeans_result, kmeans_result_prefix);
  std::vector<std::vector<size_t> > spectral_result;
  read_pairs_with_prefix(spectral_result, spectral_result_prefix);

  std::map<size_t, size_t> label_map;
  for(size_t i=0;i<spectral_result.size();++i){
    label_map.insert(std::make_pair(spectral_result[i][0], spectral_result[i][1]));
  }

  std::ofstream ofs(outfile.c_str());
  for(size_t i=0;i<kmeans_result.size();++i){
    ofs << kmeans_result[i][0] << "\t";
    ofs << label_map[kmeans_result[i][1]] << "\n";
  }

  return 0;
}

int main(int argc, char** argv) {
  std::cout << "Spectral clustering\n\n";
  time_t start, end, mid;
  std::vector<std::pair<std::string,time_t> > times;
  time(&start);

  std::string datafile;
  std::string graph_analytics_dir = "../graph_analytics/";
  std::string svd_dir = "../collaborative_filtering/";
  std::string kmeans_dir = "./";
  std::string mpi_args;
  size_t num_clusters = 0;
  size_t num_nearests = 30;
  float sigma = 1.0;
  float epsilon = 0.0;
  size_t pre_kmeans_clusters = 0;
  size_t sv = 0;
  //parse command line
  graphlab::command_line_options clopts(
          "Spectral clustering. The input data file is provided by the "
          "--data argument which is non-optional. The format of the data file is a "
          "collection of lines, where each line contains a data id followed by a "
          "comma or white-space separated list of numeric values representing a vector. "
          "Every line must have the same number of values. The required --clusters=N "
          "argument denotes the number of clusters to generate.");
  clopts.attach_option("data", datafile,
          "Input file. Each line holds a data id followed by a white-space "
          "or comma separated numeric vector");
  clopts.attach_option("clusters", num_clusters,
          "The number of clusters to create");
  clopts.attach_option("sigma", sigma,
          "Scale parameter for Gaussian kernel");
  clopts.attach_option("t-nearest", num_nearests,
          "Number of nearest neighbors (=t). Will use only the t-nearest similarities "
          "for each datapoint. If set at 0, will use all similarities.");
  clopts.attach_option("similarity-thres", epsilon,
          "Threshold to discard small similarities");
  clopts.attach_option("svd-dir", svd_dir,
          "Path to the directory where Graphlab svd is located");
  clopts.attach_option("kmeans-dir", kmeans_dir,
          "Path to the directory where Graphlab kmeans is located");
  clopts.attach_option("graph-analytics-dir", graph_analytics_dir,
          "Path to the directory where Graphlab eigen_vector_normalization is located");
  clopts.attach_option("pre-kmeans-clusters", pre_kmeans_clusters,
          "If set, will perform kmeans as a preprocess with the given cluster number.");
  clopts.attach_option("mpi-args", mpi_args,
          "If set, will execute mipexec with the given arguments. "
          "For example, --mpi-args=\"-n [N machines] --hostfile [host file]\"");
  clopts.attach_option("sv", sv,
          "Number of vectors in each iteration in the Lanczos svd.");
  if (!clopts.parse(argc, argv))
    return EXIT_FAILURE;
  if (datafile == "") {
    std::cout << "--data is not optional\n";
    return EXIT_FAILURE;
  }
  if (num_clusters == 0) {
    std::cout << "--cluster is not optional\n";
    return EXIT_FAILURE;
  }

  std::vector<std::string> remove_opts;
  remove_opts.push_back("--data");
  remove_opts.push_back("--svd-dir");
  remove_opts.push_back("--graph-analytics-dir");
  remove_opts.push_back("--kmeans-dir");
  remove_opts.push_back("--clusters");
  remove_opts.push_back("--sigma");
  remove_opts.push_back("--similarity-thres");
  remove_opts.push_back("--mpi-args");
  remove_opts.push_back("--t-nearest");
  remove_opts.push_back("--pre-kmeans-clusters");
  remove_opts.push_back("--sv");
  std::string other_args = get_arg_str_without(argc, argv, remove_opts);

  //preprocess by kmeans for fast clustering
  if(pre_kmeans_clusters > 0){
    if(pre_kmeans_clusters < num_clusters){
      std::cout << "the number of --pre-kmeans-clusters must be bigger than the number of clusters\n";
      return EXIT_FAILURE;
    }
    time(&mid);
    if(call_kmeans_as_preprocess(mpi_args, datafile, kmeans_dir, pre_kmeans_clusters, other_args) == false)
      return EXIT_FAILURE;
    //modify settings
    datafile = datafile + ".pre.centers";
    num_nearests = 0;
    time(&end);
    times.push_back(std::pair<std::string, time_t>("kmeans preprocess",(end - mid)));
  }

  //construct graph laplacian
  time(&mid);
  if (call_graph_laplacian_construction(mpi_args, datafile, sigma, epsilon,
      num_nearests, other_args) == false) {
    return EXIT_FAILURE;
  }
  time(&end);
  times.push_back(std::pair<std::string, time_t>("graph laplacian",(end - mid)));

  //eigen value decomposition
  //read number of data
  size_t num_data = 0;
  const std::string datanum_filename = datafile + ".datanum";
  std::ifstream ifs(datanum_filename.c_str());
  if (!ifs) {
    std::cout << "can't read number of data." << std::endl;
    return EXIT_FAILURE;
  }
  ifs >> num_data;
  //determine the sv of Lanczos method
  if(sv == 0){
    sv = get_lanczos_rank(num_clusters, num_data);
  }else{
    if(sv < num_clusters)
      sv = num_clusters;
  }
  time(&mid);
  if (call_svd(mpi_args, datafile, svd_dir, num_clusters, sv, num_data,
      other_args) == false) {
    return EXIT_FAILURE;
  }
  if (call_eigen_vector_normalization(mpi_args, datafile, graph_analytics_dir,
      num_clusters, sv, num_data, other_args) == false) {
    return EXIT_FAILURE;
  }
  time(&end);
  times.push_back(std::pair<std::string, time_t>("eigen decomposition",(end - mid)));

  //run kmeans
  time(&mid);
  if (call_kmeans(mpi_args, datafile, kmeans_dir, num_clusters, other_args)
      == false) {
    return EXIT_FAILURE;
  }
  time(&end);
  times.push_back(std::pair<std::string, time_t>("kmeans",(end - mid)));

  //recover cluster membership if preprocess with kmeans was done
  if(pre_kmeans_clusters > 0){
    //remove ".pre.centers"
    datafile = datafile.substr(0, datafile.size() - 12);
    recover_labels(datafile);
  }

  time(&end);

  std::cout << "computation times:\n";
  for(size_t i=0;i<times.size();++i){
    std::cout << "process " << i+1 << "\t" << times[i].first << "\t" << times[i].second << " sec\n";
  }
  std::cout << "Overall processing time of spectral clustering is " << (end - start) << " sec\n";

  return EXIT_SUCCESS;
}

