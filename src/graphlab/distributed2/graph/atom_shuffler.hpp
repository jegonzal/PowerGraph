#ifndef GRAPHLAB_ATOM_SHUFFLER_HPP
#define GRAPHLAB_ATOM_SHUFFLER_HPP



#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <omp.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/filesystem.hpp>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/stl_util.hpp>

#include <graphlab/distributed2/graph/atom_file.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {

  

  /**
   *  The atom shuffler is a distributed object that behaves like a
   *  smart DHT for atom_file construction.
   */
  template <typename VertexData, typename EdgeData>
  class atom_shuffler {
  public:
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;
    typedef atom_shuffler<VertexData, EdgeData> atom_shuffler_type;
    typedef atom_file<VertexData, EdgeData> atom_file_type;
    typedef boost::unordered_map<vertex_id_t, vertex_id_t> global2local_type;

    struct atom_info {
      atom_file_type atom_file;
      global2local_type global2local;
      // atom_info(const size_t& atomid) { 
      //   atom_file.atom_id() = atomid;
      // }
      vertex_id_t get_local_vid(vertex_id_t gvid) {
        global2local_type::iterator iter 
          = global2local.find(gvid);
        
        if(iter == global2local.end()) {
          logstream(LOG_FATAL) 
            << "Trying to find vertex " << gvid 
            << " in atom info for atom id " << atom_file.atom_id()
            << std::endl;
          ASSERT_TRUE(iter != global2local.end());
        }
        return iter->second;
      }
    };

    typedef std::vector< atom_info* > atom_map_type;    


    struct vertex_args {
      vertex_id_t vid;
      vertex_color_type vcolor;
      vertex_data_type vdata;
      vertex_args() : vid(-1), vcolor(-1) { }
      vertex_args(const vertex_id_t& vid,
                  const vertex_color_type& vcolor,
                  const vertex_data_type& vdata) :
        vid(vid), vcolor(vcolor), vdata(vdata) { } 
      void load(iarchive& iarc) {
        iarc >> vid >> vcolor >> vdata;
      }
      void save(oarchive& oarc) const {
        oarc << vid << vcolor << vdata;
      }
    }; // end of struct args

    struct edge_args {
      vertex_id_t source;
      vertex_color_type source_color;
      vertex_id_t target;
      edge_data_type edata;
      edge_args() : source(-1), source_color(-1), 
                    target(-1) { }
      edge_args(const vertex_id_t& source,
                const vertex_color_type& source_color,
                const vertex_id_t& target,
                const edge_data_type& edata) :
        source(source), source_color(source_color),
        target(target), edata(edata) { }
      void load(iarchive& iarc) {
        iarc >> source >> source_color >> target >> edata;
      }
      void save(oarchive& oarc) const {
        oarc << source << source_color << target << edata;
      }
    }; // end of struct args


    struct edge_args_bndry {
      vertex_id_t source;
      vertex_id_t target;
      vertex_color_type target_color;
      edge_args_bndry() : source(-1), target(-1) { }
      edge_args_bndry(const vertex_id_t& source,
                      const vertex_id_t& target,
                      const vertex_color_type& target_color) :
        source(source), target(target), target_color(target_color) { }
      void load(iarchive& iarc) {
        iarc >> source >>  target >> target_color;
      }
      void save(oarchive& oarc) const {
        oarc << source << target << target_color;
      }
    }; // end of struct args





  private:
    dc_dist_object< atom_shuffler_type > rmi;
    
    size_t num_atoms;
    std::vector<procid_t> vertex2atomid;

    atom_map_type atomid2info;
    std::vector<graphlab::mutex> atomid2lock;


  public:
    atom_shuffler(distributed_control& dc) : rmi(dc, this) { }



    ~atom_shuffler() {
      //      rmi.full_barrier();
      // Clear the atom info map
      for(size_t i = 0; i < atomid2info.size(); ++i) {
        atomid2lock[i].lock();
        if(atomid2info[i] != NULL) {
          delete atomid2info[i];
          atomid2info[i] = NULL;
        }
        atomid2lock[i].unlock();
      }
    } // end of destructor



      

    void add_vertex_local_vec(const std::vector<vertex_args>& args) {
      for(size_t i = 0; i < args.size(); ++i) 
        add_vertex_local(args[i].vid, args[i].vcolor, args[i].vdata);
    } // end of add vertex local vec
     


    void add_vertex_local(const vertex_id_t& gvid,
                          const vertex_color_type& vcolor,
                          const vertex_data_type& vdata) {

      const procid_t atomid(vertex2atomid.at(gvid));
      ASSERT_TRUE(is_local(atomid));
      ASSERT_LT(atomid, atomid2lock.size());
      // Get the atom info
      atom_info& ainfo(get_atom_info(atomid));
      atomid2lock[atomid].lock();
      // test to see if the vertex has been added already
      if(ainfo.global2local.find(gvid) == ainfo.global2local.end()) { 
        // first add
        const vertex_id_t localvid(ainfo.atom_file.globalvids().size());
        ainfo.atom_file.globalvids().push_back(gvid);
        ainfo.atom_file.vcolor().push_back(vcolor);
        ainfo.atom_file.atom().push_back(atomid);
        ainfo.global2local[gvid] = localvid;
        ainfo.atom_file.vdata().push_back(vdata);
      }
      atomid2lock[atomid].unlock();
    } // end of add vertex local



    void add_edge_local_vec(const std::vector<edge_args>& args) {
      for(size_t i = 0; i < args.size(); ++i) 
        add_edge_local(args[i].source, args[i].source_color,
                       args[i].target, args[i].edata);
    }// end of add edge local vec



    void add_edge_local(const vertex_id_t& source_gvid,
                        const vertex_color_type& source_vcolor,
                        const vertex_id_t& target_gvid,
                        const edge_data_type& edata) {
      ASSERT_LT(target_gvid, vertex2atomid.size());
      const procid_t atomid(vertex2atomid[target_gvid]);
      ASSERT_TRUE(is_local(atomid));
      ASSERT_LT(atomid, atomid2lock.size());
      ASSERT_NE(source_gvid, target_gvid);
      // Get the atom info
      atom_info& ainfo(get_atom_info(atomid));
      atomid2lock[atomid].lock();
      const vertex_id_t target_lvid(ainfo.get_local_vid(target_gvid)); 
      // The source may not be local yet
      if(ainfo.global2local.find(source_gvid) == ainfo.global2local.end()) {
        const vertex_id_t lvid(ainfo.atom_file.globalvids().size());
        ainfo.atom_file.globalvids().push_back(source_gvid);
        ainfo.global2local[source_gvid] = lvid;
        ainfo.atom_file.vcolor().push_back(source_vcolor);
        ainfo.atom_file.atom().push_back(vertex2atomid.at(source_gvid));
      }
      const vertex_id_t source_lvid(ainfo.get_local_vid(source_gvid));
      ASSERT_NE(source_lvid, target_lvid);
      ainfo.atom_file.edge_src_dest().push_back(std::make_pair(source_lvid, 
                                                               target_lvid));
      ainfo.atom_file.edata().push_back(edata);
      atomid2lock[atomid].unlock();
    } // end of add edge local





    void add_bndry_out_edge_vec(const std::vector<edge_args_bndry>& args) {
      for(size_t i = 0; i < args.size(); ++i) 
        add_bndry_out_edge_local(args[i].source, args[i].target, 
                                 args[i].target_color);
    }// end of add edge local vec



    void add_bndry_out_edge_local(const vertex_id_t& source_gvid,
                                  const vertex_id_t& target_gvid,
                                  const vertex_color_type& target_vcolor) {          
      ASSERT_NE(source_gvid, target_gvid);
      ASSERT_LT(source_gvid, vertex2atomid.size());
      ASSERT_LT(target_gvid, vertex2atomid.size());
      const procid_t source_atomid(vertex2atomid[source_gvid]);
      ASSERT_TRUE(is_local(source_atomid));
      const procid_t target_atomid(vertex2atomid[target_gvid]);
      ASSERT_LT(source_atomid, atomid2lock.size());
      // Get the atom info
      atom_info& ainfo(get_atom_info(source_atomid));
      ASSERT_EQ(ainfo.atom_file.atom_id(), source_atomid);
      atomid2lock[source_atomid].lock();
      const vertex_id_t source_lvid(ainfo.get_local_vid(source_gvid)); 
      // The target may not be local
      if(ainfo.global2local.find(target_gvid) == ainfo.global2local.end()) {
        const vertex_id_t lvid(ainfo.atom_file.globalvids().size());
        ainfo.atom_file.globalvids().push_back(target_gvid);
        ainfo.global2local[target_gvid] = lvid;
        ainfo.atom_file.vcolor().push_back(target_vcolor);
        ainfo.atom_file.atom().push_back(target_atomid);
      }
      const vertex_id_t target_lvid(ainfo.get_local_vid(target_gvid));
      ASSERT_NE(source_lvid, target_lvid);
      ainfo.atom_file.edge_src_dest().push_back(std::make_pair(source_lvid, 
                                                               target_lvid));
      atomid2lock[source_atomid].unlock();
    } // end of add edge local



    


    void shuffle(const std::string& atom_index_fname,
                 const std::string& partition_fname,
                 const std::string& new_atom_path,
                 const std::string& atom_prefix = "atom_") {
     
      // Load the atom index file
      atom_index_file aindex; 
      typedef std::pair<size_t, std::string> intstr_t;
      std::vector< intstr_t > local_fnames;
      { // Determine Local Filenames ==========================================
        logstream(LOG_INFO) 
          << "Loading atom index file: " << atom_index_fname << std::endl;    
        aindex.read_from_file(atom_index_fname);
        ASSERT_GT(aindex.atoms.size(), 0);
        local_fnames.resize(aindex.atoms.size());
        for(size_t i = 0; i < aindex.atoms.size(); ++i)
          local_fnames[i] = std::make_pair(i, aindex.atoms[i].file);
        // compute the fnames that are used by this machine
        std::vector< std::vector< intstr_t > > partition_fnames;
        rmi.gather_partition(local_fnames, partition_fnames);
        // update the local fnames      
        local_fnames = partition_fnames[rmi.procid()];
      }
      // Load the partitioning information
      procid_t num_atoms(-1);
      { 
        logstream(LOG_INFO) 
          << "Loading the partitioning file: " << partition_fname 
          << std::endl;    
        vertex2atomid.resize(aindex.nverts);
        std::ifstream fin(partition_fname.c_str());
        assert(fin.good());
        procid_t max_atomid(0);
        for(size_t i = 0; i < vertex2atomid.size(); ++i) {
          procid_t atomid(-1);
          fin >> atomid;
          ASSERT_TRUE(fin.good());
          ASSERT_NE(atomid, procid_t(-1));
          max_atomid = std::max(max_atomid, atomid);
          vertex2atomid[i] = atomid;
        }
        fin.close();
        num_atoms = max_atomid + 1;
        logstream(LOG_INFO) 
          << "Num atoms: " << num_atoms  << std::endl;      
      }  
              
      { // Initialize the atom files for this machine =========================
        logstream(LOG_INFO)
          << "Initializing all atom files." << std::endl;
        atomid2info.resize(num_atoms, NULL);
        atomid2lock.resize(num_atoms);
        for(procid_t atomid = 0; atomid < num_atoms; ++atomid) {         
          atomid2lock[atomid].lock();
          // if this is the owning machine setup the map and files
          if(is_local(atomid)) {
            // get the atom info (and create it)
            atom_info& ainfo(get_atom_info(atomid, true));
            std::stringstream strm;
            strm << atom_prefix
                 << std::setw(5) << std::setfill('0')
                 << atomid << ".atom";
            ainfo.atom_file.filename() = strm.str();
            ainfo.atom_file.protocol() = "file";
            // set the id of the atom_file associated with the atom info
            ainfo.atom_file.atom_id() = atomid;
          } // end of if owning machine
          atomid2lock[atomid].unlock();
        } //end of for loop
      } // end of initialize atom files for this machine
          

      rmi.full_barrier();

      std::set<vertex_id_t> allgvids;
      { // Shuffle the vertex data ============================================
        logstream(LOG_INFO) << "Shuffle vertex data." << std::endl;
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          logstream(LOG_INFO) 
            << "Loading atom file: " << local_fnames[i].second << std::endl;
          atom_file_type afile;
          afile.input_filename("file", local_fnames[i].second);
          afile.load_all();
          afile.atom_id() = local_fnames[i].first;
          //          afile.costly_checking();
          ASSERT_LE(afile.vdata().size(), afile.globalvids().size());
          ASSERT_LE(afile.vdata().size(), afile.vcolor().size());
          logstream(LOG_INFO) << "Sending vertex data." << std::endl;
          typedef vertex_args args;
          const size_t ONE_MB(size_t(1) << 20);
          const size_t BUFFER_SIZE(ONE_MB / sizeof(args));
          std::vector< std::vector<args> >  proc2buffer(rmi.numprocs());
          for(size_t j = 0, vdata_index=0; 
              j < afile.globalvids().size(); ++j) {
            const vertex_id_t gvid(afile.globalvids()[j]);
            const procid_t old_atomid( afile.atom()[j] );
            ASSERT_LT(gvid, vertex2atomid.size());
            const procid_t new_atomid( vertex2atomid[gvid] );
            ASSERT_LT(old_atomid, num_atoms);
            if(old_atomid == afile.atom_id()) {
              ASSERT_LT(j, afile.vcolor().size());
              const vertex_color_type vcolor(afile.vcolor()[j]);
              ASSERT_LT(vdata_index, afile.vdata().size());
              const vertex_data_type& vdata(afile.vdata()[vdata_index++]);
              if(is_local(new_atomid)) {
                add_vertex_local(gvid, vcolor, vdata);
              } else {
                const procid_t owner(owning_machine(new_atomid));
                proc2buffer[owner].push_back(args(gvid, vcolor, vdata));
                if(proc2buffer[owner].size() > BUFFER_SIZE) {
                  rmi.remote_call(owner, 
                                  &atom_shuffler::add_vertex_local_vec,
                                  proc2buffer[owner]);
                  proc2buffer[owner].clear();
                }
              } //end of if else
            }
          } // end of loop over vertex data
          // Flush buffers
          for(size_t i = 0; i < proc2buffer.size(); ++i) {
            if(!proc2buffer[i].empty()) {
              rmi.remote_call(i,
                              &atom_shuffler::add_vertex_local_vec,
                              proc2buffer[i]);
              proc2buffer[i].clear();
            }
          } // end of loop over flush buffers
        } // end of for loop
      } // end of shuffle vertex data

      rmi.full_barrier();

      { // Shuffle the edge data ============================================
        logstream(LOG_INFO) << "Shuffle edge data." << std::endl;
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          logstream(LOG_INFO) 
            << "Loading atom file: " << local_fnames[i].second << std::endl;
          atom_file_type afile;
          afile.input_filename("file", local_fnames[i].second);
          afile.load_all();
          afile.atom_id() = local_fnames[i].first;
          logstream(LOG_INFO) << "Sending edge data." << std::endl;
          const size_t ONE_MB(size_t(1) << 20);
          const size_t BNDRY_BUFFER_SIZE(ONE_MB / sizeof(edge_args_bndry));
          const size_t IN_BUFFER_SIZE(ONE_MB / sizeof(edge_args));
          std::vector< std::vector<edge_args> > 
            proc2buffer_in(rmi.numprocs());
          std::vector< std::vector<edge_args_bndry> > 
            proc2buffer_bndry(rmi.numprocs());
          for(size_t j = 0, edata_index = 0; 
              j < afile.edge_src_dest().size(); ++j) {
            const vertex_id_t source_lvid(afile.edge_src_dest()[j].first);
            const vertex_id_t target_lvid(afile.edge_src_dest()[j].second);
            ASSERT_NE(source_lvid, target_lvid);
            
            ASSERT_LT(source_lvid, afile.vcolor().size());
            ASSERT_LT(target_lvid, afile.vcolor().size());
            const vertex_color_type source_vcolor(afile.vcolor()[source_lvid]);
            const vertex_color_type target_vcolor(afile.vcolor()[target_lvid]);

            ASSERT_LT(source_lvid, afile.atom().size());
            ASSERT_LT(target_lvid, afile.atom().size());
            const procid_t source_old_atomid(afile.atom()[source_lvid]);
            const procid_t target_old_atomid(afile.atom()[target_lvid]);
            ASSERT_TRUE(source_old_atomid == afile.atom_id() || 
                        target_old_atomid == afile.atom_id());
            const bool HAVE_EDGE_DATA(target_old_atomid == afile.atom_id());

            ASSERT_LE(source_lvid, afile.globalvids().size());
            ASSERT_LE(target_lvid, afile.globalvids().size());
            const vertex_id_t source_gvid(afile.globalvids()[source_lvid]);
            const vertex_id_t target_gvid(afile.globalvids()[target_lvid]); 

            ASSERT_LE(target_gvid, vertex2atomid.size());
            ASSERT_LE(source_gvid, vertex2atomid.size());
            const procid_t target_atomid( vertex2atomid[target_gvid] );            
            const procid_t source_atomid( vertex2atomid[source_gvid] );        
            ASSERT_LT(source_atomid, num_atoms);
            ASSERT_LT(target_atomid, num_atoms);
            
            // If I have the edge data then I am responsible for sends
            if(HAVE_EDGE_DATA) {
              // send the edge data to the target vertex
              const edge_data_type& edata(afile.edata()[edata_index++]);
              if(is_local(target_atomid)) {
                add_edge_local(source_gvid, source_vcolor, 
                               target_gvid, edata);
              } else {
                const procid_t owner(owning_machine(target_atomid));
                proc2buffer_in[owner].push_back(edge_args(source_gvid, source_vcolor, 
                                                  target_gvid, edata));
                if(proc2buffer_in[owner].size() > IN_BUFFER_SIZE) {
                  rmi.remote_call(owner, 
                                  &atom_shuffler::add_edge_local_vec,
                                  proc2buffer_in[owner]);
                  proc2buffer_in[owner].clear();
                }
              } //end of if else
              // send the out edge to the source vertex
              if(is_local(source_atomid)) {
                add_bndry_out_edge_local(source_gvid, 
                                         target_gvid, 
                                         target_vcolor);
              } else {
                const procid_t owner(owning_machine(source_atomid));
                proc2buffer_bndry[owner].push_back(edge_args_bndry(source_gvid, 
                                                      target_gvid, 
                                                      target_vcolor));
                if(proc2buffer_bndry[owner].size() > BNDRY_BUFFER_SIZE) {
                  rmi.remote_call(owner, 
                                  &atom_shuffler::add_bndry_out_edge_vec,
                                  proc2buffer_bndry[owner]);
                  proc2buffer_bndry[owner].clear();
                }
              } //end of if else
            }
          } // end of loop over edge data
          // Flush buffers
          for(size_t i = 0; i < proc2buffer_in.size(); ++i) {
            if(!proc2buffer_in[i].empty()) {
              rmi.remote_call(i,
                              &atom_shuffler::add_edge_local_vec,
                              proc2buffer_in[i]);
            }
            proc2buffer_in[i].clear();
          } // end of loop over flush buffers
          for(size_t i = 0; i < proc2buffer_bndry.size(); ++i) {
            if(!proc2buffer_bndry[i].empty()) {
              rmi.remote_call(i,
                              &atom_shuffler::add_bndry_out_edge_vec,
                              proc2buffer_bndry[i]);
              proc2buffer_bndry[i].clear();
            }
          } // end of loop over flush buffers
        } // end of for loop
      } // end of shuffle edge data

      logstream(LOG_INFO) << "Finished shuffling!" << std::endl;

      rmi.full_barrier();

      { // Emit actual atom files =============================================
        logstream(LOG_INFO) << "Saving atom files." << std::endl;
        // Loop through each of the atoms managed locally
        for(size_t atomid = 0; atomid < atomid2info.size(); ++atomid) {
          atomid2lock[atomid].lock();
          if(atomid2info[atomid] != NULL) {
            atom_info& ainfo(get_atom_info(atomid));
            const std::string fname(new_atom_path + "/" + 
                                    ainfo.atom_file.filename());
            logstream(LOG_INFO) << "Writing final atom file "
                                << fname << std::endl;
            ainfo.atom_file.write_to_file(ainfo.atom_file.protocol(),
                                          fname);
            // Delete the pointer
            assert(atomid2info[atomid] != NULL);
            delete atomid2info[atomid];
            atomid2info[atomid] = NULL;
          } // end of if statement
          // release lock
          atomid2lock[atomid].unlock();
        } // end of for loop
      } // end of emit atom files

      rmi.full_barrier();            

      if(rmi.procid() == 0) 
        std::cout << "Finished shuffle!"
                  << std::endl;

      
    } // end of shuffle
      


  private: // helper functions
    
    procid_t owning_machine(procid_t atomid) {
      return atomid % rmi.numprocs();
    }
      
    bool is_local(procid_t atomid) {
      return owning_machine(atomid) == rmi.procid();
    }
      
    atom_info& get_atom_info(procid_t atomid, bool create = false) {
      assert(is_local(atomid));
      // Check that this is the owning machine and the file has been allocated
      assert(atomid < atomid2info.size());
      if(create) {
        assert(atomid2info[atomid] == NULL);
        atomid2info[atomid] = new atom_info();
        assert(atomid2info[atomid] != NULL);
      }
      assert(atomid2info[atomid] != NULL);
      return *atomid2info[atomid];
    }









  public:

    static void 
    rebuild_atoms_from_partitioning(const std::string& atom_index_fname,
                                    const std::string& partition_fname,
                                    const std::string& new_atom_path,
                                    const std::string& atom_prefix) {
      logstream(LOG_INFO) 
        << "Initializing distributed shuffler object." 
        << std::endl;        
      dc_init_param param;         
      if( ! init_param_from_mpi(param) ) {
        logstream(LOG_FATAL) 
          << "Failed MPI laucher!" << std::endl;
      }
      param.initstring = "buffered_send=yes, ";
      param.numhandlerthreads = 5;
      distributed_control dc(param);
      dc.full_barrier();      
      // Create the atom shuffler
      atom_shuffler atom_shuffler(dc);
      dc.full_barrier();
      atom_shuffler.shuffle(atom_index_fname, 
                            partition_fname,
                            new_atom_path,
                            atom_prefix);
      dc.full_barrier();
      logstream(LOG_INFO) << "Finished." << std::endl;      
    } // end of graph partition to atom index






  }; // End of atom shuffler





  












}; // end namespace graphlab














#include <graphlab/macros_undef.hpp>






#endif
