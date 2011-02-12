#ifndef GRAPHLAB_ATOM_SHUFFLER_HPP
#define GRAPHLAB_ATOM_SHUFFLER_HPP



#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/stl_util.hpp>

#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/distributed2/graph/partitioning/raw_fragment.hpp>



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

    struct atom_info {
      atom_file_type atom_file;
      std::map<vertex_id_t, vertex_id_t> global2local;
      std::string vdatafn;
      std::string edatafn;
      std::string atomfn;
      std::ofstream vdatastream;
      std::ofstream edatastream;        
      vertex_id_t get_local_vid(vertex_id_t gvid) {
        std::map<vertex_id_t, vertex_id_t>::iterator iter 
          = global2local.find(gvid);
        assert(iter != global2local.end());
        return iter->second;
      }
    };


    typedef std::vector< atom_info* > atom_map_type;

  private:


    dc_dist_object< atom_shuffler_type > rmi;
    size_t num_atoms;
    std::string path;
    atom_map_type atomid2info;
    std::vector<graphlab::mutex> atomid2lock;

  public:
    atom_shuffler(distributed_control& dc,
                  const size_t& num_atoms,
                  const std::string& path,
                  const std::string& atom_prefix = "atom_") : 
      rmi(dc, this), num_atoms(num_atoms), path(path),
      atomid2info(num_atoms, NULL), atomid2lock(num_atoms) { 
      // Initialize the atom files for this machine
      for(procid_t atomid = 0; atomid < num_atoms; ++atomid) {         
        atomid2lock[atomid].lock();
        // if this is the owning machine setup the map and files
        if(is_local(atomid)) {
          // get the atom info (and create it)
          atom_info& ainfo(get_atom_info(atomid, true));
          { // open vdata temporary storage file
            std::stringstream strm;
            strm << path << "/"
                 << "tmp_vdata_"
                 << std::setw(3) << std::setfill('0')
                 << atomid
                 << ".bin";
            ainfo.vdatafn = strm.str();
            ainfo.vdatastream.open(ainfo.vdatafn.c_str(), 
                                   std::ios::binary | 
                                   std::ios::out |
                                   std::ios::trunc );
            assert(ainfo.vdatastream.good());
          }

          { // open edata temporary storage file
            std::stringstream strm;
            strm << path << "/"
                 << "tmp_edata_"
                 << std::setw(3) << std::setfill('0')
                 << atomid
                 << ".bin";
            ainfo.edatafn = strm.str();
            ainfo.edatastream.open(ainfo.edatafn.c_str(), 
                                   std::ios::binary | 
                                   std::ios::out |
                                   std::ios::trunc );
            assert(ainfo.edatastream.good());
          }

          { // determine atom filename
            std::stringstream strm;
            strm << path << "/"
                 << atom_prefix
                 << std::setw(3) << std::setfill('0')
                 << atomid
                 << ".atom";
            ainfo.atomfn = strm.str();
          }
          // set the id of the atom_file associated with the atom info
          ainfo.atom_file.atom_id() = atomid;
        } //end of if owning machine
        atomid2lock[atomid].unlock();
      } //end of for loop
      dc.services().full_barrier();
    } //end of constructor
      


    ~atom_shuffler() {
      rmi.services().full_barrier();
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
      

    void add_vertex_local(const procid_t& to_atomid,
                          const vertex_id_t& gvid,
                          const procid_t& atomid,
                          const vertex_color_type& vcolor,
                          const vertex_data_type& vdata) {
      assert(is_local(to_atomid));
      atomid2lock[to_atomid].lock();
      // Get the atom info
      atom_info& ainfo(get_atom_info(to_atomid));
      // test to see if the vertex has been added already
      if(ainfo.global2local.find(gvid) == ainfo.global2local.end()) { // first add
        vertex_id_t localvid = ainfo.atom_file.globalvids().size();
        ainfo.atom_file.globalvids().push_back(gvid);
        ainfo.atom_file.vcolor().push_back(vcolor);
        ainfo.atom_file.atom().push_back(atomid);
        ainfo.global2local[gvid] = localvid;
        assert(ainfo.vdatastream.good());
        {
          oarchive oarc(ainfo.vdatastream);
          oarc << vdata;
        }
      }
      atomid2lock[to_atomid].unlock();
    }

    void add_vertex(const procid_t& to_atomid,
                    const vertex_id_t& gvid,
                    const procid_t& atomid,
                    const vertex_color_type& vcolor,
                    const vertex_data_type& vdata) {

      if(is_local(to_atomid)) {
        // local add
        add_vertex_local(to_atomid, gvid, atomid, vcolor, vdata);
      } else {
        // remote add
        rmi.remote_call(owning_machine(to_atomid),
                        &atom_shuffler_type::add_vertex_local,
                        to_atomid, gvid, atomid, vcolor, vdata);
      }
    } // end of add_vertex



    void add_edge_local(const procid_t& to_atomid,
                        const vertex_id_t& source_gvid,
                        const vertex_id_t& target_gvid,
                        const edge_data_type& edata) {
      assert(is_local(to_atomid));
      atomid2lock[to_atomid].lock();
      // Get the atom info
      atom_info& ainfo(get_atom_info(to_atomid));
      vertex_id_t source_lvid(ainfo.get_local_vid(source_gvid));
      vertex_id_t target_lvid(ainfo.get_local_vid(target_gvid));        
      ainfo.atom_file.edge_src_dest().push_back(std::make_pair(source_lvid, 
                                                               target_lvid));
      assert(ainfo.edatastream.good());
      {  // THIS SCOPE IS ABSOLUTELY NECESSARY:
        // to ensure that the oarchive does not flush outside the
        // mutex.
        oarchive oarc(ainfo.edatastream);
        oarc << edata;      
      }
      assert(ainfo.edatastream.good());
      atomid2lock[to_atomid].unlock();
    } // end of add edge local


    void add_edge(const procid_t& to_atomid,
                  const vertex_id_t& source_gvid,
                  const vertex_id_t& target_gvid,
                  const edge_data_type& edata) {
      if(is_local(to_atomid)) {
        add_edge_local(to_atomid, source_gvid, target_gvid, edata);
      } else {
        rmi.remote_call(owning_machine(to_atomid),
                        &atom_shuffler_type::add_edge_local,
                        to_atomid, source_gvid, target_gvid, edata);

      }
    } // end of add edge



    void load_vertex_data(const std::vector< std::string >& fnames,
                          const std::vector< procid_t >& vertex2atomid,
                          const std::vector< vertex_color_type >& vertex2color) {
      foreach(const std::string fname, fnames) {
        raw_fragment frag;
        { // Get the local structure description
          std::string absfname = path + "/" + fname;
          std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);
          iarc >> frag;
          fin.close();
        }
        { // read all the vdata
          std::string absfname = path + "/" + frag.vertex_data_filename;
          std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);        
          // Loop through all the structures reading in the vertex data
          for(vertex_id_t i = 0, vid = frag.begin_vertex; 
              vid < frag.end_vertex; ++vid, ++i) {
            assert(vid < vertex2atomid.size());
            assert(vid < vertex2color.size());
            vertex_data_type vdata;
            assert(fin.good());
            iarc >> vdata;
            assert(fin.good());
            // send the vertex data to the owning atoms
            add_vertex(vertex2atomid[vid],
                       vid,
                       vertex2atomid[vid],
                       vertex2color[vid],
                       vdata);
            // Loop over the vertex neighbors and send to neighbor atoms
            foreach(vertex_id_t neighbor_vid, frag.neighbor_ids[i]) {
              assert(neighbor_vid < vertex2atomid.size());
              // if the neighbor is stored in a different atom file then
              // send the vertex data to the neighbor for ghosting purposes
              if(vertex2atomid[neighbor_vid] != vertex2atomid[vid]) {
                // send the vertex data to the owning atoms
                add_vertex(vertex2atomid[neighbor_vid],
                           vid,
                           vertex2atomid[vid],
                           vertex2color[vid],
                           vdata);
              } // end of if neighbor is in different atom
            }
          }
          fin.close();
        } // end of read all vdata
      } // end of adding all vertices


    } // end of load vertex data



    void load_edge_data(const std::vector< std::string >& fnames,
                        const std::vector< procid_t >& vertex2atomid) {    
      foreach(const std::string fname, fnames) {
        raw_fragment frag;    
        { // Get the local structure description
          std::string absfname = path + "/" + fname;
          std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);
          iarc >> frag;
          fin.close();
        }

        { // read all the edge data
          std::string absfname = path + "/" + frag.edge_data_filename;
          std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);        
          // Loop through all the structures reading in the vertex data
          for(vertex_id_t i = 0, target = frag.begin_vertex; 
              target < frag.end_vertex; ++target, ++i) {
            assert(target < vertex2atomid.size());
            assert(i < frag.in_neighbor_ids.size());
            foreach(vertex_id_t source, frag.in_neighbor_ids[i]) {
              assert(source < vertex2atomid.size());
              edge_data_type edata;
              assert(fin.good());
              iarc >> edata;
              assert(fin.good());
              add_edge(vertex2atomid[source],
                       source, target, edata);            
              if(vertex2atomid[source] != vertex2atomid[target])
                add_edge(vertex2atomid[target],
                         source, target, edata);                       
            } // end of loop over out neighbors
          } // end of loop over vertices in file
          fin.close();
        } // end of read all vdata
      } // end of for loop


    } // end of load edge data






    void emit_atoms() {
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
        }
        // release lock
        atomid2lock[atomid].unlock();
      } // end of for loop
    } // end of emit atoms
    

  private: // helper functions
    
    procid_t owning_machine(procid_t atomid) {
      return atomid % rmi.numprocs();
    }
      
    bool is_local(procid_t atomid) {
      return owning_machine(atomid) == rmi.procid();
    }
      
    atom_info& get_atom_info(procid_t atomid, bool create = false) {
      // Check that this is the owning machine and the file has been allocated
      assert(atomid < atomid2info.size());
      if(create) {
        assert(atomid2info[atomid] == NULL);
        atomid2info[atomid] = new atom_info();
      }
      assert(atomid2info[atomid] != NULL);
      return *atomid2info[atomid];
    }









  public:

    static void build_atoms_from_partitioning(const std::string& path) {      
      std::cout << "Initializing distributed communication layer" << std::endl;
      dc_init_param param;      
      if ( !init_param_from_env(param) ) {
        std::cout << "Failed to get environment variables" << std::endl;
        assert(false);
      }
      distributed_control dc(param);
      dc.services().full_barrier();
    
      if(dc.procid() == 0) 
        std::cout << "Loading vertex coloring and partitioning." 
                  << std::endl;

      // load the maps
      vertex_color_type max_color = 0;
      std::vector<vertex_color_type> vertex2color;
      {
        std::string absfname = path + "/coloring.txt";
        std::ifstream fin(absfname.c_str());
        while(fin.good()) {
          vertex_color_type vcolor;
          fin >> vcolor;
          max_color = std::max(max_color, vcolor);
          if(fin.good()) vertex2color.push_back(vcolor);
        }
        fin.close();
      }
      size_t num_colors( max_color + 1 );


      if(dc.procid() == 0) 
        std::cout << "Coloring: " << num_colors
                  << std::endl;


      procid_t max_atomid = 0;
      std::vector<procid_t> vertex2atomid;
      {
        std::string absfname = path + "/partitioning.txt";
        std::ifstream fin(absfname.c_str());
        while(fin.good()) {
          procid_t atomid;
          fin >> atomid;
          max_atomid = std::max(max_atomid, atomid);
          if(fin.good()) vertex2atomid.push_back(atomid);
        }
        fin.close();
      }
      size_t num_atoms( max_atomid + 1 );
      if(dc.procid() == 0) 
        std::cout << "Num atoms: " << num_atoms
                  << std::endl;



      // get the filenames accessible by this process
      std::vector<std::string> fnames;
      {
        raw_fragment::list_structure_files(path, fnames);
        // compute the fnames that are used by this machine
        std::vector< std::vector< std::string > > partition_fnames;
        dc.services().gather_partition(fnames, partition_fnames);

        // update the local fnames      
        fnames = partition_fnames[dc.procid()];
        std::cout << "Assigned Names: " << std::endl;
        foreach(std::string fname, fnames) 
          std::cout << "(" << fname << ")" << '\t';
        std::cout << std::endl;
      }

      
      if(dc.procid() == 0) { 
        std::cout << "Initializing distributed shuffler object.  ";
        std::cout.flush();
      }
      // Create the atom shuffler
      atom_shuffler atom_shuffler(dc, 
                                  num_atoms, 
                                  path);

      dc.services().full_barrier();
      if(dc.procid() == 0) 
        std::cout << "Finished." << std::endl;

      {
        if(dc.procid() == 0) {
          std::cout << "Loading all vertex data from graph fragments. ";
          std::cout.flush();
        }
        atom_shuffler.load_vertex_data(fnames, 
                                       vertex2atomid, 
                                       vertex2color );
        dc.services().full_barrier();
        if(dc.procid() == 0) 
          std::cout << "Finished." << std::endl;
      }
    
      {
        if(dc.procid() == 0)  {
          std::cout << "Loading all edge data from graph fragments. ";        
          std::cout.flush();
        }

        atom_shuffler.load_edge_data(fnames, vertex2atomid);
        dc.services().full_barrier();   
        if(dc.procid() == 0)  
          std::cout << "Finished." << std::endl;
      }

      {
        if(dc.procid() == 0) {
          std::cout << "Emitting all atoms. ";
          std::cout.flush();
        }
        // Build the actual atom files
        atom_shuffler.emit_atoms();
        dc.services().full_barrier();   
        
        if(dc.procid() == 0) 
          std::cout << "Finished." << std::endl;
      }
        

   
    } // end of graph partition to atom index






  }; // End of atom shuffler












































}; // end namespace graphlab
#include <graphlab/macros_undef.hpp>






#endif
