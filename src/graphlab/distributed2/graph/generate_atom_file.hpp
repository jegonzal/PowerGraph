/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GENERATE_ATOM_FILE_HPP
#define GENERATE_ATOM_FILE_HPP
#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {
  
  template<typename VertexData, typename EdgeData,
           typename VFunction, typename EFunction>
  void generate_atom_files(const std::string& path,
                           const std::string& atompath,
                           VFunction vfun, EFunction efun) {
    typedef atom_file<VertexData, EdgeData> atom_file_type;
    const std::string vdata_suffix(".vdata_txt");
    const std::string edata_suffix(".edata_txt");
    dc_init_param param;         
    if( ! init_param_from_mpi(param) ) {
      logstream(LOG_FATAL) 
        << "Failed MPI laucher!" << std::endl;
    }
    distributed_control dc(param);
    dc.full_barrier();      
    logstream(LOG_INFO) 
      << "Initializing distributed shuffler object." 
      << std::endl;    

    // Determine the local vdata files on this machine
    std::vector<std::string> local_fnames;
    fs_util::list_files_with_suffix(path, vdata_suffix,
                                    local_fnames);
    std::vector< std::vector<std::string> > partition_fnames;
    dc.gather_partition(local_fnames, partition_fnames);
    local_fnames = partition_fnames.at(dc.procid());   
    std::vector< procid_t > local_atomids(local_fnames.size(), 
                                          procid_t(-1));

    std::vector<procid_t> atomid2proc;
    { // comute local atom ids
      procid_t starting_atomid(0);
      for(size_t i = 0; i < dc.procid(); ++i)
        starting_atomid += partition_fnames.at(i).size();
      for(size_t i = 0; i < local_fnames.size(); ++i)
        local_atomids.at(i) = starting_atomid + i;      
      for(size_t i = 0; i < partition_fnames.size(); ++i) {
        for(size_t j = 0; j < partition_fnames[i].size(); ++j) {
          atomid2proc.push_back(i);
        }
      }
    }

    std::vector< std::string > local_atom_fnames(local_fnames.size());
    { // compute atom file names
      for(size_t i = 0; i < local_atom_fnames.size(); ++i) 
        local_atom_fnames.at(i) =  
          fs_util::change_suffix(local_fnames[i], ".atom");        
    }

    typedef std::pair<vertex_id_t, vertex_id_t> edge_t;
    typedef std::pair<procid_t, edge_t> atom_edge_t;
    std::vector< std::vector<atom_edge_t> > proc2out_edges(dc.numprocs());

    typedef std::map<procid_t, std::vector<vertex_id_t> > 
      atomid2vertexids_type;
    std::vector< atomid2vertexids_type > gather_vec(dc.numprocs());    
    atomid2vertexids_type& atomid2vertexids(gather_vec.at(dc.procid()));
    // Build out each atom file
    for(size_t i = 0; i < local_fnames.size(); ++i) {
      const std::string vdata_fname = path + "/" +
        fs_util::change_suffix(local_fnames[i], vdata_suffix);
      const std::string edata_fname = path + "/" +
        fs_util::change_suffix(local_fnames[i], edata_suffix);
      atom_file_type afile;
      afile.load_from_txt(vdata_fname, edata_fname, vfun, efun);
      afile.atom_id() = local_atomids[i];
      { // update the atomid2vertexids map
        std::vector<vertex_id_t>& 
          localverts(atomid2vertexids[afile.atom_id()]);
        localverts.resize(afile.globalvids().size());
        for(size_t j = 0; j < localverts.size(); ++j)
          localverts[j] = afile.globalvids().at(j);
        foreach(const edge_t& edge, afile.edge_src_dest()) {
          const procid_t atomid(afile.atom()[edge.first]);
          const vertex_id_t source(afile.globalvids()[edge.first]);
          const vertex_id_t target(afile.globalvids()[edge.second]);
          const edge_t global_edge(std::make_pair(source, target));
          if(atomid != afile.atom_id()) 
            proc2out_edges[atomid2proc[atomid]].
              push_back(std::make_pair(atomid, global_edge));
        }

        
      }
      afile.filename() = local_atom_fnames.at(i);
      afile.protocol() = "file";
      // construct the atom filename
      const std::string fname(path + "/" + afile.filename() +"_tmp");
      logstream(LOG_INFO) << "Saving atom file: " 
                          << fname
                          << std::endl;
      afile.write_to_file("file", fname );
    }

    std::vector< std::vector<edge_t> > atom2out_edges(atomid2proc.size());
    {
      std::vector< std::vector<atom_edge_t> > result;
      mpi_tools::all2all(proc2out_edges, result);
      for(size_t i = 0; i < result.size(); ++i) {
        for(size_t j = 0; j < result[i].size(); ++j) {
          atom2out_edges[result[i][j].first].push_back(result[i][j].second);
        }
      }
    }


    const size_t ROOT_NODE(0);
    dc.gather(gather_vec, ROOT_NODE);
    std::vector<procid_t> vid2atomid;
    if(dc.procid() == 0) {
      logstream(LOG_INFO) << "Computing vertex to atom map."  << std::endl;
      // compute max vertex id
      vertex_id_t max_vid(0);
      typedef atomid2vertexids_type::value_type pair_type;
      for(size_t i = 0; i < gather_vec.size(); ++i) {         
        foreach(const pair_type& pair, gather_vec[i]) {
          for(size_t j = 0; j < pair.second.size(); ++j) {
            max_vid = std::max(max_vid, pair.second[j]);
          }
        }
      }
      const vertex_id_t nverts(max_vid + 1);
      logstream(LOG_INFO) << "Num vertices " << nverts << std::endl;
      // Construct inverse map
      vid2atomid.resize(nverts, procid_t(-1));
      // Invert the map
      for(size_t i = 0; i < gather_vec.size(); ++i) {
        foreach(const pair_type& pair, gather_vec[i]) {
          const procid_t atomid(pair.first);
          for(size_t j = 0; j < pair.second.size(); ++j) {
            const vertex_id_t vid(pair.second[j]);
            ASSERT_LT(vid, nverts);
            ASSERT_EQ(vid2atomid[vid], vertex_id_t(-1));
            vid2atomid[vid] = atomid;
          }
        }
      }
      // check all the entries
      for(size_t i = 0; i < vid2atomid.size(); ++i)
        ASSERT_NE(vid2atomid[i], procid_t(-1));
      // broadcast the map to all machines
      dc.broadcast(vid2atomid, true);
    } else {
      dc.broadcast(vid2atomid, false);
    }

    // Loop through local atom files
    logstream(LOG_INFO) << "Updating local atoms."  << std::endl;
    



    // Build out each atom file
    for(size_t i = 0; i < local_atom_fnames.size(); ++i) {
      const std::string tmp_fname(path + "/" + local_atom_fnames[i] + "_tmp");
      logstream(LOG_INFO) << "Reloading " << tmp_fname << std::endl;
      atom_file_type afile;
      afile.input_filename("file", tmp_fname);
      afile.load_all();
      std::map<vertex_id_t, vertex_id_t> global2local;
      // Update the atom location information
      for(size_t j = 0; j < afile.globalvids().size(); ++j) {
        const vertex_id_t gvid(afile.globalvids()[j]);
        global2local[gvid] = j;
        ASSERT_LT(gvid, vid2atomid.size());
        const procid_t atomid(vid2atomid[gvid]);
        ASSERT_NE(atomid, procid_t(-1));
        afile.atom()[j] = atomid;
      }
      // Adde the extra edges
      for(size_t j = 0; j < atom2out_edges[i].size(); ++j) {
        const edge_t edge(atom2out_edges[i][j]);
        ASSERT_TRUE(global2local.find(edge.first) != global2local.end());
        const vertex_id_t source_lvid(global2local[edge.first]);        
        if(global2local.find(edge.second) == global2local.end()) {
          const vertex_id_t target_lvid = afile.globalvids().size();
          global2local[edge.second] = target_lvid;
          afile.globalvids().push_back(edge.second);          
        }
        const vertex_id_t target_lvid(global2local[edge.second]);
        afile.edge_src_dest().push_back
          (std::make_pair(source_lvid, target_lvid));

      }

      const std::string fname( atompath + "/" + local_atom_fnames[i]);
      // Resave the atom file
      logstream(LOG_INFO) << "Final save of atom file: " 
                          << fname
                          << std::endl;
      afile.write_to_file("file", fname );
    }

    

  } // end of generate atom file

}
#include <graphlab/macros_undef.hpp>
#endif
