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

bool call_graph_laplacian(const std::string& mpi_args,
    const std::string& filename, const std::string& format,
    const bool normalized_cut, const bool ratio_cut, const std::string& args) {
  std::stringstream strm;
  if(mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << "./graph_laplacian ";
  strm << " --graph=" << filename;
  strm << " --format=" << format;
//  strm << " --normalized-cut=" << normalized_cut;
//  strm << " --ratio-cut=" << ratio_cut;
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
  if(mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << svd_dir << "svd " + filename + ".glap";
  strm << " --rows=" << num_data+1;
  strm << " --cols=" << num_data;
  strm << " --nsv=" << num_clusters;
  strm << " --nv=" << rank;
  strm << " --max_iter=4";
  strm << " --quiet=1";
  strm << " --save_vectors=1";
  strm << " --ortho_repeats=3";
  //strm << " --id=1";
  //strm << " --prediction=" << filename + ".";
  strm << " --prediction=" << filename;
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
    const std::string& filename, const size_t num_clusters, const size_t rank,
    const size_t num_data, const std::string& args) {
  std::stringstream strm;
  if(mpi_args.length() > 0)
    strm << "mpiexec " << mpi_args << " ";
  strm << "./eigen_vector_normalization";
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
  if(mpi_args.length() > 0)
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

//select good rank
int get_lanczos_rank(const size_t num_clusters, const size_t num_data) {
  size_t rank = 1;
  if (num_data < 1000) {
    if (num_clusters + 10 <= num_data)
      rank = num_clusters + 10;
    else
      rank = num_data;
  } else if (num_data < 10000) {
    rank = num_clusters + 25;
  } else if (num_data < 100000) {
    rank = num_clusters + 50;
  } else if (num_data < 1000000) {
    rank = num_clusters + 80;
  } else {
    rank = num_clusters + 100;
  }
  return rank;
//  return num_clusters + 1;
}

int main(int argc, char** argv) {
  std::cout << "Graph partitioning (normalized cut)\n\n";

  std::string graph_dir;
  std::string format = "adj";
  std::string svd_dir = "../collaborative_filtering/";
  std::string kmeans_dir = "../clustering/";
  std::string mpi_args;
  size_t num_partitions = 2;
  bool normalized_cut = true;
  bool ratio_cut = false;
  size_t sv = 0;
  //parse command line
  graphlab::command_line_options clopts(
          "Graph partitioning (normalized cut)");
  clopts.attach_option("graph", graph_dir,
                       "The graph file. This is not optional. Vertex ids must start from 1 "
                       "and must not skip any numbers.");
  clopts.attach_option("format", format,
                       "The graph file format. If \"weight\" is set, the program will read "
                       "the data file where each line holds [id1] [id2] [weight].");
  clopts.attach_option("partitions", num_partitions,
                       "The number of partitions to create");
  clopts.attach_option("svd-dir", svd_dir,
                       "Path to the directory of Graphlab svd");
  clopts.attach_option("kmeans-dir", kmeans_dir,
                       "Path to the directory of Graphlab kmeans");
  clopts.attach_option("mpi-args", mpi_args,
                       "If set, will execute mipexec with the given arguments. "
                       "For example, --mpi-args=\"-n [N machines] --hostfile [host file]\"");
  clopts.attach_option("sv", sv,
                       "Number of vectors in each iteration in the Lanczos svd.");
//  clopts.attach_option("normalized-cut", normalized_cut,
//                       "do normalized cut");
//  clopts.attach_option("ratio-cut", ratio_cut,
//                       "do ratio cut");
  if (!clopts.parse(argc, argv))
    return EXIT_FAILURE;
  if (graph_dir == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }
//  if(normalized_cut == true && ratio_cut == true){
//    std::cout << "Both normalized-cut and ratio-cut are true. Ratio cut is selected.\n";
//    normalized_cut = false;
//  }else if(normalized_cut == false && ratio_cut == false){
//    std::cout << "Both normalized-cut and ratio-cut are false. Ratio cut is selected.\n";
//    ratio_cut = true;
//  }
  std::vector<std::string> remove_opts;
  remove_opts.push_back("--graph");
  remove_opts.push_back("--format");
  remove_opts.push_back("--svd-dir");
  remove_opts.push_back("--kmeans-dir");
  remove_opts.push_back("--partitions");
  remove_opts.push_back("--mpi-args");
  remove_opts.push_back("--sv");
//  remove_opts.push_back("--normalized-cut");
//  remove_opts.push_back("--ratio-cut");
  std::string other_args = get_arg_str_without(argc, argv, remove_opts);

  //construct graph laplacian
  if (call_graph_laplacian(mpi_args, graph_dir, format, normalized_cut,
      ratio_cut, other_args) == false) {
    return EXIT_FAILURE;
  }

  //eigen value decomposition
  //read number of data
  size_t num_data = 0;
  const std::string datanum_filename = graph_dir + ".datanum";
  std::ifstream ifs(datanum_filename.c_str());
  if (!ifs) {
    std::cout << "can't read number of data." << std::endl;
    return false;
  }
  ifs >> num_data;
  //determine the rank of Lanczos method
  if(sv == 0){
    sv = get_lanczos_rank(num_partitions, num_data);
  }else{
    if(sv < num_partitions)
      sv = num_partitions;
  }
  if (call_svd(mpi_args, graph_dir, svd_dir, num_partitions, sv, num_data,
      other_args) == false) {
    return EXIT_FAILURE;
  }
  if (call_eigen_vector_normalization(mpi_args, graph_dir, num_partitions, sv,
      num_data, other_args) == false) {
    return EXIT_FAILURE;
  }

  //kmeans
  if (call_kmeans(mpi_args, graph_dir, kmeans_dir, num_partitions, other_args)
      == false) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

