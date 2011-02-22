#ifndef GRAPHLAB_ATOM_SHUFFLER_HPP
#define GRAPHLAB_ATOM_SHUFFLER_HPP



#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include <graphlab/util/timer.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/stl_util.hpp>

#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>



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

    size_t num_colors;
    std::vector<vertex_color_type> vertex2color;
    
    size_t num_atoms;
    std::vector<procid_t> vertex2atomid;

    
    adjacency_list alist;
    std::vector< vertex_id_t > vertex2proc;
    std::vector< std::set< procid_t > >  neighbor_atoms;
    std::vector< graphlab::mutex > neighbor_locks;

    atom_map_type atomid2info;
    std::vector<graphlab::mutex> atomid2lock;



  public:
    atom_shuffler(distributed_control& dc) : rmi(dc, this) { }



    ~atom_shuffler() {
      rmi.full_barrier();
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



    void add_neighbor_atom(const vertex_id_t vid, 
                           const procid_t neighbor_atom) {
      // if vid is stored locally
      if(vertex2proc[vid] == rmi.procid() ) {
        add_neighbor_atom_local(vid, neighbor_atom);
      } else {
        // remote add
        rmi.remote_call(vertex2proc[vid],
                        &atom_shuffler_type::add_neighbor_atom_local,
                        vid, neighbor_atom);
      }
    } // end of add atom neighbor


    void add_neighbor_atom_local(const vertex_id_t vid, 
                                 const procid_t neighbor_atom) {
      assert(vid < vertex2proc.size());
      assert(vertex2proc[vid] == rmi.procid());
      assert(neighbor_atom < num_atoms);
      // convert to the local vertex id address
      assert(alist.global2local.find(vid) != alist.global2local.end());
      vertex_id_t localvid = alist.global2local[vid];

      assert(localvid < neighbor_atoms.size());
      assert(localvid < neighbor_locks.size());
      neighbor_locks[localvid].lock();
      neighbor_atoms[localvid].insert(neighbor_atom);
      neighbor_locks[localvid].unlock();
    } // end of add atom neighbor local

      

    void add_vertex_local(const procid_t& to_atomid,
                          const vertex_id_t& gvid,
                          const procid_t& atomid,
                          const vertex_color_type& vcolor,
                          const vertex_data_type& vdata) {
      assert(is_local(to_atomid));
      assert(to_atomid < atomid2lock.size());
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
        { // for scoping of archive object
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
      { 
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

    


    void shuffle(const std::string& path,
                 const std::string& atom_path,
                 const std::string& atom_prefix = "atom_") {
     
      

      std::vector<std::string> local_fnames;
      { // Determine Local Filenames ==========================================
        if(rmi.procid() == 0) 
          std::cout << "Computing local filenames."
                    << std::endl;    
        adjacency_list::list_vlist_files(path, local_fnames);
        // compute the fnames that are used by this machine
        std::vector< std::vector< std::string > > partition_fnames;
        rmi.gather_partition(local_fnames, partition_fnames);
        // update the local fnames      
        local_fnames = partition_fnames[rmi.procid()];
     
      }
  
      { // Load local adjacency lists fragments ===============================
        if(rmi.procid() == 0) 
          std::cout << "Loading local graph structures."
                    << std::endl;
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          std::string absfname = path + "/" + local_fnames[i];
          // std::cout << "(" <<  rmi.procid() << " - " 
          //           << local_fnames[i] << ")" << std::endl;
          alist.load(absfname);
        }
      }

      size_t nverts(0), nedges(0);
      {  // compute the total number of vertices and edges
        std::vector<size_t> gather_count(rmi.numprocs(), 0);
        gather_count.at(rmi.procid()) = alist.local_vertices.size();
        rmi.all_gather(gather_count);
        for(size_t i = 0; i < gather_count.size(); ++i) 
          nverts += gather_count.at(i);
        for(size_t i = 0; i < alist.in_neighbor_ids.size(); ++i) 
          nedges += alist.in_neighbor_ids.at(i).size();
        gather_count.at(rmi.procid()) = nedges;
        rmi.all_gather(gather_count);
        nedges = 0;
        for(size_t i = 0; i < gather_count.size(); ++i) 
          nedges += gather_count.at(i);
        std::cout << "Nvertices: " << nverts << std::endl
                  << "Nedges:    " << nedges << std::endl;
      }





      { // Load the partitioning information =========================
        if(rmi.procid() == 0) 
          std::cout << "Loading partitioning."
                    << std::endl;      
        size_t max_atomid(0);
        std::string absfname = path + "/partitioning.txt";
        std::ifstream fin(absfname.c_str());
        if(!fin.good()) { // randomly partition if not possible
          std::cout << "Using random partitioning." << std::endl;
          vertex2atomid.resize(nverts, procid_t(-1));
          for(size_t i = 0; i < vertex2atomid.size(); ++i) { 
            size_t atomid = (i % 1000);
            max_atomid = std::max(max_atomid, atomid);            
            vertex2atomid[i] = atomid;
          }
        } else {
          vertex2atomid.reserve(nverts);
          while(fin.good()) {
            size_t atomid;
            fin >> atomid;
            if(fin.good()) {
              assert(atomid < std::numeric_limits<procid_t>::max());
              max_atomid = std::max(max_atomid, atomid);            
              vertex2atomid.push_back(procid_t(atomid));
            }
          }
        }
        fin.close();
        num_atoms = max_atomid + 1;        
        if(rmi.procid() == 0) {
          std::cout << "Num atoms: " << num_atoms
                    << std::endl;      
          // Save the partitioning in the shared atom folder
          absfname = atom_path + "/partitioning.txt";
          std::ofstream fout(absfname.c_str());
          assert(fout.good());
          for(size_t i = 0; i < vertex2atomid.size(); ++i) {
            fout << size_t(vertex2atomid[i]) << '\n';
            assert(fout.good());
          }
          fout.close();
        }
      }

      { // Load the coloring information =============================
        if(rmi.procid() == 0) 
          std::cout << "Loading vertex coloring."
                    << std::endl;      
        size_t max_color(0);
        std::string absfname = path + "/coloring.txt";
        std::ifstream fin(absfname.c_str());
        if(!fin.good()) {// file does not exist 
          vertex2color.resize(vertex2atomid.size(), 0);
        }
        while(fin.good()) {
          size_t vcolor;
          fin >> vcolor;
          if(fin.good()) {
            assert(vcolor < 
                   std::numeric_limits<vertex_color_type>::max());
            max_color = std::max(max_color, vcolor);
            vertex2color.push_back(vertex_color_type(vcolor));
          }
        }
        fin.close();
        num_colors = max_color + 1;
        if(rmi.procid() == 0)  {
          std::cout << "Coloring: " << num_colors
                    << std::endl;
          
          // Save the coloring in the shared atom folder
          absfname = atom_path + "/coloring.txt";
          std::ofstream fout(absfname.c_str());
          assert(fout.good());
          for(size_t i = 0; i < vertex2color.size(); ++i) {
            fout << size_t(vertex2color[i]) << '\n';
            assert(fout.good());
          }
          fout.close();
        }
      }




      {
        // check the vertices in the alists
        for(size_t i = 0; i < alist.local_vertices.size(); ++i) {
          assert(alist.local_vertices[i] < vertex2atomid.size());
        }
        // resize auxiliarary datastructures
        assert(alist.in_neighbor_ids.size() == alist.local_vertices.size());
        assert(alist.global2local.size() == alist.local_vertices.size());
        neighbor_atoms.resize(alist.local_vertices.size());
        neighbor_locks.resize(alist.local_vertices.size());
      }


      assert(vertex2atomid.size() == vertex2color.size());
      if(rmi.procid() == 0) 
        std::cout << "Vertices: " << vertex2color.size() << std::endl;
      



      timer ti;
      ti.start();

      rmi.full_barrier();

      { // Compute the vertex2proc map ======================================
        if(rmi.procid() == 0) 
          std::cout << "Computing vertex to proc map"
                    << std::endl;
        const size_t ROOT_NODE(0);
        std::vector< std::vector< vertex_id_t > > proc2vertices(rmi.numprocs());
        proc2vertices[rmi.procid()] = alist.local_vertices;
        rmi.gather(proc2vertices, ROOT_NODE);
        // if this is the root node then invert the map
        if(rmi.procid() == ROOT_NODE) {
          vertex2proc.resize(vertex2atomid.size(), procid_t(-1));
          for(size_t i = 0; i < proc2vertices.size(); ++i) {
            for(size_t j = 0; j < proc2vertices[i].size(); ++j) { 
              vertex_id_t vid = proc2vertices[i][j];
              if(vid >= vertex2atomid.size()) 
                std::cout << "Error!: " << vid << std::endl;
              assert(vid < vertex2atomid.size());
              assert(vertex2proc[vid] == procid_t(-1));
              vertex2proc[vid] = i;
            }
          }
          // send the complete map
          rmi.broadcast(vertex2proc, true);
        } else {
          // receive the complete map
          rmi.broadcast(vertex2proc, false);
        }
        assert(vertex2proc.size() == vertex2atomid.size());
        // check the final map
        for(size_t i = 0; i < vertex2proc.size(); ++i) {
          assert(vertex2proc[i] < rmi.numprocs());
        }
      } // end of compute initial vertex2proc map


      
      rmi.full_barrier();
      
      { // compute neighbor atoms for all local vertices ======================
        if(rmi.procid() == 0) 
          std::cout << "Computing atom neighbors for each vertex."
                    << std::endl;
        for(size_t i = 0; i < alist.in_neighbor_ids.size(); ++i) {
          assert(i < alist.local_vertices.size());
          vertex_id_t target = alist.local_vertices[i];
          assert(target < vertex2atomid.size());
          for(size_t j = 0; j < alist.in_neighbor_ids[i].size(); ++j) {
            vertex_id_t source = alist.in_neighbor_ids[i][j];
            assert(source < vertex2atomid.size());
            add_neighbor_atom(target, vertex2atomid[source]);
            add_neighbor_atom(source, vertex2atomid[target]);
          }
        }
      } // end of compute neighbor atoms for all vertices


      rmi.full_barrier();
        
      { // Initialize the atom files for this machine =========================
        if(rmi.procid() == 0) 
          std::cout << "Initializing all atom files."
                    << std::endl;
        atomid2info.resize(num_atoms, NULL);
        atomid2lock.resize(num_atoms);
        for(procid_t atomid = 0; atomid < num_atoms; ++atomid) {         
          atomid2lock[atomid].lock();
          // if this is the owning machine setup the map and files
          if(is_local(atomid)) {
            // get the atom info (and create it)
            atom_info& ainfo(get_atom_info(atomid, true));
            { // open vdata temporary storage file
              std::stringstream strm;
              strm << path << "/" << "tmp_vdata_"
                   << std::setw(5) << std::setfill('0')
                   << atomid << ".bin";
              ainfo.vdatafn = strm.str();
              ainfo.vdatastream.open(ainfo.vdatafn.c_str(), 
                                     std::ios::binary | std::ios::out |
                                     std::ios::trunc );
              assert(ainfo.vdatastream.good());
            }
            
            { // open edata temporary storage file
              std::stringstream strm;
              strm << path << "/" << "tmp_edata_"
                 << std::setw(5) << std::setfill('0')
                 << atomid << ".bin";
              ainfo.edatafn = strm.str();
              ainfo.edatastream.open(ainfo.edatafn.c_str(), 
                                     std::ios::binary | std::ios::out |
                                     std::ios::trunc );
              assert(ainfo.edatastream.good());
            }

            { // determine atom filename
              std::stringstream strm;
              strm << atom_path << "/"  << atom_prefix
                 << std::setw(5) << std::setfill('0')
                 << atomid << ".atom";
              ainfo.atomfn = strm.str();
            }
            // set the id of the atom_file associated with the atom info
            ainfo.atom_file.atom_id() = atomid;
          } // end of if owning machine
          atomid2lock[atomid].unlock();
        } //end of for loop
      } // end of initialize atom files for this machine
          

      rmi.full_barrier();


      { // Shuffle the vertex data ============================================
        if(rmi.procid() == 0) 
          std::cout << "Shuffling vertex data."
                    << std::endl;
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
            assert(fin.good());
            assert(localvid < alist.local_vertices.size());
            vertex_id_t vid = alist.local_vertices[localvid];
            assert(vid < vertex2atomid.size());
            assert(vid < vertex2color.size());
            add_vertex(vertex2atomid[vid],
                       vid,
                       vertex2atomid[vid],
                       vertex2color[vid],
                       vdata);
            assert(localvid < neighbor_atoms.size());
            foreach(procid_t neighbor_atom, neighbor_atoms[localvid]) {
              assert(neighbor_atom < num_atoms);
              // if the neighbor is stored in a different atom file
              // then send the vertex data to the neighbor for
              // ghosting purposes
              if(neighbor_atom != vertex2atomid[vid]) {
                // send the vertex data to the owning atoms
                add_vertex(neighbor_atom,
                           vid,
                           vertex2atomid[vid],
                           vertex2color[vid],
                           vdata);
              } // end of if neighbor is in different atom
            } // end of loop over neighbors
            localvid++; // successful add so increment the local vid counter
          } // end of loop over single vertex data file
          fin.close();
        } // end of loop over all vertex data files
        assert(localvid == alist.local_vertices.size());        
      } // end of shuffle the vertex data


      rmi.comm_barrier();
      rmi.barrier();

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
            assert(localvid < alist.in_neighbor_ids.size());
            vertex_id_t target(alist.local_vertices[localvid]);
            assert(target < vertex2atomid.size());
            // try to read in all the neighbors
            for(size_t j = 0; j < alist.in_neighbor_ids[localvid].size(); ++j) {
              vertex_id_t source(alist.in_neighbor_ids[localvid][j]);
              assert(source < vertex2atomid.size());
              edge_data_type edata;
              iarc >> edata;
              assert(fin.good());
              add_edge(vertex2atomid[source], source, target, edata);
              if(vertex2atomid[source] != vertex2atomid[target])
                add_edge(vertex2atomid[target], source, target, edata); 
            }
            localvid++;
            fin.peek();
          } // end of while loop
          fin.close();
        } // end of for loop
        assert(localvid == alist.local_vertices.size());
      } // end of shuffle the edge data

      std::cout << "Entering final barrier on " << rmi.procid() << std::endl;
      rmi.comm_barrier();
      rmi.barrier();
      rmi.dc().full_barrier();
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
      std::cout << "Initializing distributed communication layer" << std::endl;
      dc_init_param param;      
      if ( !init_param_from_env(param) ) {
        std::cout << "Failed to get environment variables" << std::endl;
        assert(false);
      }
      param.initstring = "buffered_send=yes";
      param.numhandlerthreads = 30;
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
