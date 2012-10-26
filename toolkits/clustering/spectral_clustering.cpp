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
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <graphlab.hpp>

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
    const std::string& args) {
  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << "./graph_laplacian_for_sc ";
  strm << " --data=" << filename;
  strm << " --sigma=" << sigma;
  strm << " --similarity-thres=" << epsilon;
  strm << " " << args;
  std::cout << "CALLING >" << strm.str() << std::endl;
  int sys_ret = system(strm.str().c_str());
  if (sys_ret != 0) {
    std::cout << "system call fails" << std::endl;
    return false;
  }

  return true;
}

bool call_svd(const std::string& mpi_args, const std::string& filename,
    const std::string& svd_dir, const size_t num_clusters, const size_t rank,
    const size_t num_data, const std::string& args) {

  std::stringstream strm;
  if (mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << svd_dir << "svd " + filename + ".glap";
  strm << " --rows=" << num_data;
  strm << " --cols=" << num_data;
  strm << " --nsv=" << rank;
  strm << " --nv=" << rank;
  strm << " --max_iter=4";
  strm << " --quiet=1";
  strm << " --save_vectors=1";
  strm << " --ortho_repeats=3";
  strm << " --id=1";
  strm << " --prediction=" << filename << ".";
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

//select good value of rank (TODO)
int get_lanczos_rank(const size_t num_clusters, const size_t num_data) {
  size_t rank = 1;
  if (num_data < 1000) {
    if (num_clusters + 10 <= num_data)
      rank = num_clusters + 40;
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
}

int main(int argc, char** argv) {
  std::cout << "Spectral clustering\n\n";

  std::string datafile;
  std::string graph_analytics_dir = "../graph_analytics/";
  std::string svd_dir = "../collaborative_filtering/";
  std::string kmeans_dir = "./";
  std::string mpi_args;
  size_t num_clusters = 0;
  float sigma = 0.1;
  float epsilon = 0.0;
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
  clopts.attach_option("graph-analytics-dir", graph_analytics_dir,
          "Path to the directory of Graphlab graph analytics tools");
  clopts.attach_option("svd-dir", svd_dir,
          "Path to the directory of Graphlab svd");
  clopts.attach_option("kmeans-dir", kmeans_dir,
          "Path to the directory of Graphlab kmeans");
  clopts.attach_option("sigma", sigma,
          "Scale parameter for Gaussian kernel");
  clopts.attach_option("similarity-thres", epsilon,
          "Threshold to discard small similarities");
  clopts.attach_option("mpi-args", mpi_args,
          "If set, will execute mipexec with the given arguments. "
          "For example, --mpi-args=\"-n [N machines] --hostfile [host file]\"");
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
  std::string other_args = get_arg_str_without(argc, argv, remove_opts);

  //construct graph laplacian
  if (call_graph_laplacian_construction(mpi_args, datafile, sigma, epsilon,
      other_args) == false) {
    return EXIT_FAILURE;
  }

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
  //determine the rank of Lanczos method
  size_t rank = get_lanczos_rank(num_clusters, num_data);
  if (call_svd(mpi_args, datafile, svd_dir, num_clusters, rank, num_data,
      other_args) == false) {
    return EXIT_FAILURE;
  }
  if (call_eigen_vector_normalization(mpi_args, datafile, graph_analytics_dir,
      num_clusters, rank, num_data, other_args) == false) {
    return EXIT_FAILURE;
  }

  //run kmeans
  if (call_kmeans(mpi_args, datafile, kmeans_dir, num_clusters, other_args)
      == false) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

