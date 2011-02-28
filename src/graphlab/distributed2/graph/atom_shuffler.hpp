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
      vertex_id_t get_local_vid(vertex_id_t gvid) {
        global2local_type::iterator iter 
          = global2local.find(gvid);
        assert(iter != global2local.end());
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
      vertex_id_t target;
      edge_data_type edata;
      edge_args() : source(-1), target(-1) { }
      edge_args(const vertex_id_t& source,
                const vertex_id_t& target,
                const edge_data_type& edata) :
        source(source), target(target), edata(edata) { }
      void load(iarchive& iarc) {
        iarc >> source >> target >> edata;
      }
      void save(oarchive& oarc) const {
        oarc << source << target << edata;
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
        info.vdata().push_back(vdata);
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
      const procid_t atomid(vertex2atomid.at(target_gvid));
      ASSERT_TRUE(is_local(atomid));
      ASSERT_LT(atomid, atomid2lock.size());
      ASSERT_NE(source_gvid, target_gvid);
      // Get the atom info
      atom_info& ainfo(get_atom_info(atomid));
      atomid2lock[atomid].lock();
      const vertex_id_t target_lvid(ainfo.get_local_vid(target_gvid)); 
      // The source may not be local yet
      if(ainfo.global2local.find(source_gvid) == ainfo.global2local.end()) {
        const vertex_id_t lvid(ainfo.atom_file.globalids().size());
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
      atomid2lock[to_atomid].unlock();
    } // end of add edge local



    


    void shuffle(const std::string& atom_index_fname,
                 const std::string& partition_fname,
                 const std::string& new_atom_path,
                 const std::string& atom_prefix = "atom_") {
     
      // Load the atom index file
      atom_index_file aindex; 
      std::vector<std::string> local_fnames;
      { // Determine Local Filenames ==========================================
        logstream(LOG_INFO) 
          << "Loading atom index file: " << atom_index_fname << std::endl;    
        aindex.read_from_file(atom_index_fname);
        ASSERT_GT(aindex.atoms.size(), 0);
        local_fnames.resize(aindex.atoms.size());
        for(size_t i = 0; i < aindex.atoms.size(); ++i)
          local_fnames[i] = aindex.atoms[i].file;
        // compute the fnames that are used by this machine
        std::vector< std::vector< std::string > > partition_fnames;
        rmi.gather_partition(local_fnames, partition_fnames);
        // update the local fnames      
        local_fnames = partition_fnames[rmi.procid()];
      }
      // Load the partitioning information
      procid_t num_atoms(-1);
      { 
        logstream(LOG_INFO) 
          << "Loading the partitioning file: " << partitioning_fname 
          << std::endl;    
        vertex2atomid.resize(index.nverts);
        std::ifstream fin(partitioning_fname.c_str());
        assert(fin.good());
        procid_t max_atomid(0);
        for(size_t i = 0; i < vertex2atomid.size(); ++i) {
          procid_t atomid(-1);
          fin >> atomid;
          ASSERT(fin.good());
          ASSERT_NE(atomid, procid_t(-1));
          max_atomid = std::max(max_atomid, atomid);
          vertex2atomid[i] = atomid;
        }
        fin.close();
        num_atoms = max_atom + 1;
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

      /// TODO: FINISH SHUFFLE CODE HERE
      { // Shuffle the vertex data ============================================
        if(rmi.procid() == 0) 
          std::cout << "Shuffling vertex data."
                    << std::endl;
        typedef vertex_args args;
        const size_t ONE_MB(size_t(1) << 20);
        const size_t BUFFER_SIZE(ONE_MB / sizeof(args));
        std::vector< std::vector<args> >  proc2buffer(rmi.numprocs());

        size_t localvid(0);
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          // get the vertex data filename from the structure filename
          std::string 
            vdata_fname(adjacency_list::
                        change_suffix(local_fnames[i],
                                      adjacency_list::vdata_suffix));
          std::string absfname = path + '/' + vdata_fname;
          std::cout << "(" <<  rmi.procid() << " - " 
                    << absfname << ")" << std::endl;

          std::ifstream fin(absfname.c_str(), 
                            std::ios::binary | std::ios::in);
          graphlab::iarchive iarc(fin);
          vertex_data_type vdata;
          for(iarc >> vdata; fin.good(); iarc >> vdata) {
            ASSERT_LT(localvid, alist.local_vertices.size());
            const vertex_id_t vid( alist.local_vertices[localvid] );
            ASSERT_LT(vid, vertex2atomid.size());
            const procid_t atomid( vertex2atomid[vid] );
            if(is_local(atomid)) {
              add_vertex_local(vid, vdata);
            } else {
              const procid_t owner(owning_machine(atomid));
              proc2buffer[owner].push_back(args(vid, vdata));
              if(proc2buffer[owner].size() > BUFFER_SIZE) {
                rmi.remote_call(owner, 
                                &atom_shuffler::add_vertex_local_vec,
                                proc2buffer[owner]);
                proc2buffer[owner].clear();
              }
            }
            localvid++; // successful add so increment the local vid counter
          } // end of loop over single vertex data file
          fin.close();
        } // end of loop over all vertex data files
        assert(localvid == alist.local_vertices.size());        
        // Flush buffers
        for(size_t i = 0; i < proc2buffer.size(); ++i) {
          if(!proc2buffer[i].empty()) {
            rmi.remote_call(i,
                            &atom_shuffler::add_vertex_local_vec,
                            proc2buffer[i]);
          }
        }
      } // end of shuffle the vertex data


      rmi.full_barrier();


      { // Shuffle the edge data ==============================================
        if(rmi.procid() == 0) 
          std::cout << "Shuffling edge data."
                    << std::endl;
        
        size_t localvid(0);
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          // get the vertex data filename from the structure filename
          std::string 
            edata_fname(adjacency_list::
                        change_suffix(local_fnames[i],
                                      adjacency_list::edata_suffix));
          std::string absfname = path + '/' + edata_fname;
          std::cout << "(" <<  rmi.procid() << " - " 
                    << absfname << ")" << std::endl;

          std::ifstream fin(absfname.c_str(), 
                            std::ios::binary | std::ios::in);
          graphlab::iarchive iarc(fin);
          assert(fin.good());
          fin.peek();
          while(fin.good()) {
            assert(localvid < alist.in_nbr_ids.size());
            vertex_id_t target(alist.local_vertices[localvid]);
            assert(target < vertex2atomid.size());
            // try to read in all the nbrs
            for(size_t j = 0; j < alist.in_nbr_ids[localvid].size(); ++j) {
              vertex_id_t source(alist.in_nbr_ids[localvid][j]);
              assert(source < vertex2atomid.size());
              edge_data_type edata;
              iarc >> edata;
              assert(fin.good());
              add_edge(vertex2atomid[source], source, target, edata);
              // if(vertex2atomid[source] != vertex2atomid[target])
              //   add_edge(vertex2atomid[target], source, target, edata); 
            }
            localvid++;
            fin.peek();
          } // end of while loop
          fin.close();
        } // end of for loop
        assert(localvid == alist.local_vertices.size());
      } // end of shuffle the edge data

      std::cout << "Entering final barrier on " << rmi.procid() << std::endl;

      rmi.full_barrier();
      std::cout << "Leaving final barrier on " << rmi.procid() << std::endl;


      if(rmi.procid() == 0)
        std::cout << "Finished shuffling in " << ti.current_time() 
                  << " seconds." << std::endl;


      { // Emit actual atom files =============================================
        if(rmi.procid() == 0) 
          std::cout << "Saving atom files."
                    << std::endl;
        namespace fs = boost::filesystem;
        // Loop through each of the atoms managed locally
        for(size_t atomid = 0; atomid < atomid2info.size(); ++atomid) {
          atomid2lock[atomid].lock();
          if(atomid2info[atomid] != NULL) {
            atom_info& ainfo(*atomid2info[atomid]);
            { // read all vdata into atom file
              ainfo.vdatastream.flush();
              ainfo.vdatastream.close();
              assert(!ainfo.vdatastream.is_open());
              std::ifstream fin(ainfo.vdatafn.c_str(), 
                                std::ios::binary | 
                                std::ios::in);
              assert(fin.good());
              iarchive iarc(fin);
              atom_file_type& afile(ainfo.atom_file);
              afile.vdata().resize(afile.globalvids().size());
              for(size_t i = 0; i < afile.vdata().size(); ++i) {
                assert(fin.good());
                iarc >> afile.vdata()[i];                                
              }
              fin.close();
              fs::path path(ainfo.vdatafn);
              fs::remove(path);
            } // end of read all vdata
            { // read all edata into atom file
              ainfo.edatastream.flush();
              ainfo.edatastream.close();
              assert(!ainfo.edatastream.is_open());
              std::ifstream fin(ainfo.edatafn.c_str(), 
                                std::ios::binary | 
                                std::ios::in);
              assert(fin.good());
              iarchive iarc(fin);
              atom_file_type& afile(ainfo.atom_file);
              afile.edata().resize(afile.edge_src_dest().size());
              for(size_t i = 0; i < afile.edata().size(); ++i) {
                assert(fin.good());
                edge_data_type edata;
                iarc >> edata;
                afile.edata()[i] = edata;                                
              }
              fin.close();
              fs::path path(ainfo.edatafn);
              fs::remove(path);
            } // end of read all edata
            
            // Save the atom
            ainfo.atom_file.write_to_file("file", ainfo.atomfn);
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
    build_atoms_from_partitioning(const std::string& path,
                                  const std::string& atom_path) {      
      std::cout << "Initializing distributed communication layer "
                << "using environment variables."  << std::endl;
      dc_init_param param;      
      // This works with BOTH mpi and rpcexec
      // if ( !init_param_from_env(param) ) {
      //   std::cout << "Failed to get environment variables." << std::endl;
      //   std::cout << "Trying MPI launcher." << std::endl;
      // } else 

      if( ! init_param_from_mpi(param) ) {
        std::cout << "Failed MPI laucher!" << std::endl;
        assert(false);
      }
      param.initstring = "buffered_send=yes, ";
      param.numhandlerthreads = 5;
      distributed_control dc(param);
      dc.full_barrier();      
      if(dc.procid() == 0) { 
        std::cout << "Initializing distributed shuffler object.  ";
        std::cout.flush();
      }
      // Create the atom shuffler
      atom_shuffler atom_shuffler(dc);
      dc.full_barrier();
      atom_shuffler.shuffle(path, atom_path);
      dc.full_barrier();
      if(dc.procid() == 0) 
        std::cout << "Finished." << std::endl;      
    } // end of graph partition to atom index






  }; // End of atom shuffler





  












}; // end namespace graphlab














#include <graphlab/macros_undef.hpp>






#endif
