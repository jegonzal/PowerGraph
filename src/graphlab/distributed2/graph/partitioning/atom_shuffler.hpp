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

    struct fraginfo {
      raw_fragment frag;
      std::vector< std::set< procid_t > >  neighbor_atoms;
      std::vector< graphlab::mutex > locks;
    };

    typedef std::map<vertex_id_t, fraginfo> id2fraginfo_type;



  private:


    dc_dist_object< atom_shuffler_type > rmi;

    size_t num_colors;
    std::vector<vertex_color_type> vertex2color;
    
    size_t num_atoms;
    std::vector<procid_t> vertex2atomid;

    std::vector<std::string> local_fnames;
    id2fraginfo_type id2fraginfo;
    std::map<vertex_id_t, vertex_id_t> fragid2proc;

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



    void add_neighbor_atom(const vertex_id_t fragid,
                           const vertex_id_t vid, 
                           const procid_t neighbor_atom) {
      assert(fragid2proc.find(fragid) != fragid2proc.end());
      procid_t vid_proc = fragid2proc[fragid];
      assert(vid_proc < rmi.numprocs());
      // if v1 is stored locally add the v2 atom to v1
      if(vid_proc == rmi.procid()) {
        add_neighbor_atom_local(fragid, vid, neighbor_atom);
      } else {
        // remote add
        rmi.remote_call(vid_proc,
                        &atom_shuffler_type::add_neighbor_atom_local,
                        fragid, vid, neighbor_atom);
      }
    } // end of add atom neighbor


    void add_neighbor_atom_local(const vertex_id_t fragid,
                                 const vertex_id_t vid, 
                                 const procid_t neighbor_atom) {
      assert(fragid2proc.find(fragid) != fragid2proc.end());
      assert(fragid2proc[fragid] == rmi.procid());
      assert(id2fraginfo.find(fragid) != id2fraginfo.end());
      assert(neighbor_atom < num_atoms);
      fraginfo& finfo(id2fraginfo[fragid]);
      assert(finfo.frag.begin_vertex <= vid);
      assert(vid < finfo.frag.end_vertex);
      vertex_id_t localvid = finfo.frag.local_id(vid);
      assert(localvid < finfo.frag.num_local_verts);
      finfo.locks.at(localvid).lock();
      finfo.neighbor_atoms.at(localvid).insert(neighbor_atom);
      finfo.locks.at(localvid).unlock();
    } // end of add atom neighbor local

      

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

    


    void shuffle(const std::string& path,
                 const std::string& atom_path,
                 const std::string& atom_prefix = "atom_") {
      { // Load the coloring information =============================
        if(rmi.procid() == 0) 
          std::cout << "Loading vertex coloring."
                    << std::endl;      
        size_t max_color(0);
        std::string absfname = path + "/coloring.txt";
        std::ifstream fin(absfname.c_str());
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



      { // Load the partitioning information =========================
        if(rmi.procid() == 0) 
          std::cout << "Loading partitioning."
                    << std::endl;      
        size_t max_atomid(0);
        std::string absfname = path + "/partitioning.txt";
        std::ifstream fin(absfname.c_str());
        while(fin.good()) {
          size_t atomid;
          fin >> atomid;
          if(fin.good()) {
            assert(atomid < std::numeric_limits<procid_t>::max());
            max_atomid = std::max(max_atomid, atomid);            
            vertex2atomid.push_back(procid_t(atomid));
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

 
      { // Determine Local Filenames ==========================================
        if(rmi.procid() == 0) 
          std::cout << "Computing local filenames."
                    << std::endl;    
        raw_fragment::list_structure_files(path, local_fnames);
        // compute the fnames that are used by this machine
        std::vector< std::vector< std::string > > partition_fnames;
        rmi.gather_partition(local_fnames, partition_fnames);
        // update the local fnames      
        local_fnames = partition_fnames[rmi.procid()];
      }
  
      { // Load local raw fragments ===========================================
        if(rmi.procid() == 0) 
          std::cout << "Loading local graph structures."
                    << std::endl;
        for(size_t i = 0; i < local_fnames.size(); ++i) {
          std::string absfname = path + "/" + local_fnames[i];
          std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);
          raw_fragment frag;
          iarc >> frag;
          id2fraginfo[frag.id].frag = frag;
          fin.close();
          id2fraginfo[frag.id].neighbor_atoms.resize(frag.num_local_verts);
          id2fraginfo[frag.id].locks.resize(frag.num_local_verts);
        }      
      }

      timer ti;
      ti.start();

      { // compute fragid2proc map ==============================================
        if(rmi.procid() == 0) 
          std::cout << "Computing frag to processor map."
                    << std::endl;
        typedef typename id2fraginfo_type::value_type id2fraginfo_pair_t;
        foreach(const id2fraginfo_pair_t& pair, id2fraginfo) 
          fragid2proc[pair.first] = rmi.procid();
        std::vector< std::map<vertex_id_t, vertex_id_t> > 
          all_fragid2proc(rmi.numprocs());
        all_fragid2proc[rmi.procid()] = fragid2proc;
        rmi.all_gather(all_fragid2proc);
        for(size_t i = 0; i < all_fragid2proc.size(); ++i) {
          if(i != rmi.procid()) {
            fragid2proc.insert(all_fragid2proc[i].begin(), 
                               all_fragid2proc[i].end());
          }
        }
      }    
      
      rmi.full_barrier();
      
      { // compute neighbor atoms for all frags ===============================
        if(rmi.procid() == 0) 
          std::cout << "Computing atom neighbors for each vertex."
                    << std::endl;
        typedef typename id2fraginfo_type::value_type id2fraginfo_pair_t;
        foreach(const id2fraginfo_pair_t& pair, id2fraginfo) {
          const raw_fragment& frag(pair.second.frag); 
          for(vertex_id_t i = 0, vid = frag.begin_vertex; 
              vid < frag.end_vertex; ++vid, ++i) {
            foreach(vertex_id_t neighbor_vid, frag.in_neighbor_ids.at(i)) {
              vertex_id_t neighbor_frag_id = 
                frag.owning_fragment(neighbor_vid);
              add_neighbor_atom(frag.id, 
                                vid, 
                                vertex2atomid.at(neighbor_vid));
              add_neighbor_atom(neighbor_frag_id, 
                                neighbor_vid, 
                                vertex2atomid.at(vid));
            }
          }
        }
      } // end of compute neighbor atoms for all frags


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
                   << std::setw(3) << std::setfill('0')
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
                 << std::setw(3) << std::setfill('0')
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
                 << std::setw(3) << std::setfill('0')
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
        typedef typename id2fraginfo_type::value_type id2fraginfo_pair_t;
        foreach(const id2fraginfo_pair_t& pair, id2fraginfo) {
          const fraginfo& finfo(pair.second);
          std::string absfname = path + "/" + 
            finfo.frag.vertex_data_filename;
          std::ifstream fin(absfname.c_str(), 
                            std::ios::binary | std::ios::in);
          assert(fin.good());
          graphlab::iarchive iarc(fin);        
          // Loop through all the structures reading in the vertex data
          for(vertex_id_t i = 0, vid = finfo.frag.begin_vertex; 
              vid < finfo.frag.end_vertex; ++vid, ++i) {
            assert(vid < vertex2atomid.size());
            assert(vid < vertex2color.size());
            vertex_data_type vdata;
            assert(fin.good());
            iarc >> vdata;
            assert(fin.good());
            // send the vertex data to the owning atoms
            add_vertex(vertex2atomid[vid],
                       vid,
                       vertex2atomid.at(vid),
                       vertex2color.at(vid),
                       vdata);
            // Loop over the vertex neighbors and send to neighbor atoms
            foreach(procid_t neighbor_atom, finfo.neighbor_atoms[i]) {
              assert(neighbor_atom < num_atoms);
              // if the neighbor is stored in a different atom file then
              // send the vertex data to the neighbor for ghosting purposes
              if(neighbor_atom != vertex2atomid[vid]) {
                // send the vertex data to the owning atoms
                add_vertex(neighbor_atom,
                           vid,
                           vertex2atomid.at(vid),
                           vertex2color.at(vid),
                           vdata);
              } // end of if neighbor is in different atom
            } // and of loop over neighbor atoms
          } // end of loop over vertices in this frag
          fin.close();
        } // end of loop over all raw fragments on this machine
      } // end of shuffle the vertex data
      rmi.comm_barrier();
      rmi.barrier();

      rmi.full_barrier();


      { // Shuffle the edge data ==============================================
        if(rmi.procid() == 0) 
          std::cout << "Shuffling edge data."
                    << std::endl;
        typedef typename id2fraginfo_type::value_type id2fraginfo_pair_t;
        foreach(const id2fraginfo_pair_t& pair, id2fraginfo) {
          const raw_fragment& frag(pair.second.frag);       
          std::string absfname = path + "/" + frag.edge_data_filename;
          std::ifstream fin(absfname.c_str(), 
                            std::ios::binary | std::ios::in);
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
        } // end of for loop
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
