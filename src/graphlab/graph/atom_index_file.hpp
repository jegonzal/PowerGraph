/**  
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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE
#define GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE

#include <string>
#include <vector>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>

namespace graphlab {


  struct atom_file_descriptor{
    
    std::string protocol;
    std::string file;
    size_t nverts;
    size_t nedges;
    std::vector<size_t> adjatoms;
    std::vector<size_t> optional_weight_to_adjatoms;
    void save(oarchive &oarc) const{
      oarc << protocol
           << file
           << nverts
           << nedges
           << adjatoms
           << optional_weight_to_adjatoms;
    }
    
    void load(iarchive &iarc) {
      iarc >> protocol
           >> file
           >> nverts
           >> nedges
           >> adjatoms
           >> optional_weight_to_adjatoms;
    }
  };

  /**
     The atom index file is a file describing the locations and
     adjacency structure of a collection of atom files. It can be
     generated from a list of atom files using the disk_graph. The
     atom index file is stored as text for human readability.
   */
  struct atom_index_file {

    size_t nverts, nedges, ncolors, natoms;
    std::vector<atom_file_descriptor> atoms;

    void read_from_file(std::string indexfile);
    void write_to_file(std::string outfilename);
  };


  /**
     Assign the atoms to parts (machines) by striping over the nparts.
     That is atom $i$ is assigned to machine 
     \f$ ret[i] = i \mod \text{nparts} \f$
   */
  std::vector<std::vector<size_t> > 
  partition_atoms_sliced(const atom_index_file& atomindex, size_t nparts);

  /**
     This function assigns atoms to each processor by partitioning the
     atom graph using METIS.
   */
  std::vector<std::vector<size_t> >
  partition_atoms(const atom_index_file& atomindex, size_t nparts);


} // end namespace graphlab

#endif

