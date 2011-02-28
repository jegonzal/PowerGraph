#ifndef GRAPHLAB_ALIST_SHUFFLER_HPP
#define GRAPHLAB_ALIST_SHUFFLER_HPP



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
#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {

  

  /**
   *  The atom shuffler is a distributed object that behaves like a
   *  smart DHT for atom_file construction.
   */
  template <typename VertexData, typename EdgeData>
  class alist_shuffler {
  public:
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;
    typedef alist_shuffler<VertexData, EdgeData> alist_shuffler_type;
    typedef atom_file<VertexData, EdgeData> atom_file_type;
    typedef boost::unordered_map<vertex_id_t, vertex_id_t> global2local_type;

    struct atom_info {
      atom_file_type atom_file;
      global2local_type global2local;
      std::string vdatafn;
      std::string edatafn;
      std::string atomfn;
      std::ofstream vdatastream;
      std::ofstream edatastream;        
      vertex_id_t get_local_vid(vertex_id_t gvid) {
        global2local_type::iterator iter 
          = global2local.find(gvid);
        assert(iter != global2local.end());
        return iter->second;
      }
    };

    typedef std::vector< atom_info* > atom_map_type;    

    struct nbr_args {
      vertex_id_t vid;
      procid_t nbr_atom;
      nbr_args(const vertex_id_t& vid = 0, 
               const procid_t& nbr_atom = 0) :
        vid(vid), nbr_atom(nbr_atom) { }
      void load(iarchive& iarc) {
        iarc >> vid >> nbr_atom;
      }
      void save(oarchive& oarc) const {
        oarc << vid << nbr_atom;
      }
    }; // end of struct args

    struct vertex_args {
      vertex_id_t vid;
      //      procid_t atomid;
      //      vertex_color_type vcolor;
      vertex_data_type vdata;
      vertex_args() : vid(-1) { }
      vertex_args(const vertex_id_t vid,
                  const vertex_data_type& vdata) :
        vid(vid), vdata(vdata) { } 
      void load(iarchive& iarc) {
        iarc >> vid >> vdata;
      }
      void save(oarchive& oarc) const {
        oarc << vid << vdata;
      }
    }; // end of struct args





  private:
    dc_dist_object< alist_shuffler_type > rmi;

    size_t num_colors;
    std::vector<vertex_color_type> vertex2color;
    
    size_t num_atoms;
    std::vector<procid_t> vertex2atomid;

    
    adjacency_list alist;
    std::vector< vertex_id_t > vertex2proc;
    std::vector< boost::unordered_set< procid_t > >  nbr_atoms;
    std::vector< graphlab::mutex > nbr_locks;

    atom_map_type atomid2info;
    std::vector<graphlab::mutex> atomid2lock;




  public:
    alist_shuffler(distributed_control& dc) : rmi(dc, this) { }



    ~alist_shuffler() {
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





 
    
    
    void add_nbr_atom_local_vec(const std::vector<nbr_args>& args) {
      for(size_t i = 0; i < args.size(); ++i)
        add_nbr_atom_local(args[i].vid, args[i].nbr_atom);
    } // end of add atom nbr local


    void add_nbr_atom_local(const vertex_id_t vid, 
                            const procid_t nbr_atom) {
      ASSERT_LT(vid, vertex2proc.size());
      ASSERT_EQ(vertex2proc[vid], rmi.procid());
      ASSERT_LT(nbr_atom, num_atoms);
      // convert to the local vertex id address
      assert(alist.global2local.find(vid) !=
             alist.global2local.end());
      const vertex_id_t localvid( alist.global2local[vid] );

      ASSERT_LT(localvid, nbr_atoms.size());
      ASSERT_LT(localvid, nbr_locks.size());
      nbr_locks[localvid].lock();
      nbr_atoms[localvid].insert(nbr_atom);
      nbr_locks[localvid].unlock();
    } // end of add atom nbr local

      

    void add_vertex_local_vec(const std::vector<vertex_args>& args) {
      for(size_t i = 0; i < args.size(); ++i) 
        add_vertex_local(args[i].vid, args[i].vdata);
    }
     


    void add_vertex_local(const vertex_id_t& gvid,
                          const vertex_data_type& vdata) {
      const procid_t atomid(vertex2atomid.at(gvid));
      const vertex_color_type vcolor(vertex2color.at(gvid));
      ASSERT_TRUE(is_local(atomid));
      ASSERT_LT(atomid, atomid2lock.size());
      atomid2lock[atomid].lock();
      // Get the atom info
      atom_info& ainfo(get_atom_info(atomid));
      // test to see if the vertex has been added already
      if(ainfo.global2local.find(gvid) == ainfo.global2local.end()) { 
        // first add
        const vertex_id_t localvid(ainfo.atom_file.globalvids().size());
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
      atomid2lock[atomid].unlock();
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
                        &alist_shuffler_type::add_vertex_local,
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
                        &alist_shuffler_type::add_edge_local,
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
          std::cout << "(" <<  rmi.procid() << " - " 
                    << local_fnames[i] << ")" << std::endl;
          alist.load(absfname);
        }
        std::cout << "CPU " << rmi.procid() << " loaded "
                  << alist.local_vertices.size() << " vertices."
                  << std::endl;
      }

      size_t nverts(0), nedges(0);
      {  // compute the total number of vertices and edges
        std::vector<size_t> gather_count(rmi.numprocs(), 0);
        gather_count.at(rmi.procid()) = alist.local_vertices.size();
        rmi.all_gather(gather_count);
        for(size_t i = 0; i < gather_count.size(); ++i) 
          nverts += gather_count.at(i);
        for(size_t i = 0; i < alist.in_nbr_ids.size(); ++i) 
          nedges += alist.in_nbr_ids.at(i).size();
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
        std::cout << "Checking all local values."
                  << std::endl;
        // check the vertices in the alists
        for(size_t i = 0; i < alist.local_vertices.size(); ++i) {
          ASSERT_LT(alist.local_vertices[i], vertex2atomid.size());
        }
        // resize auxiliarary datastructures
        assert(alist.in_nbr_ids.size() == alist.local_vertices.size());
        assert(alist.global2local.size() == alist.local_vertices.size());
        nbr_atoms.resize(alist.local_vertices.size());
        nbr_locks.resize(alist.local_vertices.size());
        assert(vertex2atomid.size() == vertex2color.size());
      }
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
      
      { // compute nbr atoms for all local vertices ======================
        if(rmi.procid() == 0) 
          std::cout << "Computing atom nbrs for each vertex."
                    << std::endl;

        typedef nbr_args args;
        const size_t nthreads(8);
        const size_t ONE_MB(size_t(1) << 20);
        const size_t BUFFER_SIZE(ONE_MB / sizeof(args));
        std::vector< std::vector<args> > 
          proc2buffer(nthreads * rmi.numprocs());
#pragma omp parallel for num_threads(nthreads) 
        for(long i = 0; i < long(alist.in_nbr_ids.size()); ++i) {
          const int threadid = omp_get_thread_num();
          const vertex_id_t target(alist.local_vertices.at(i));
          const procid_t target_proc(vertex2proc.at(target));
          const procid_t target_atom(vertex2atomid.at(target));
          const size_t target_idx(target_proc + threadid * rmi.numprocs());
          for(size_t j = 0; j < alist.in_nbr_ids[i].size(); ++j) {
            const vertex_id_t source(alist.in_nbr_ids[i][j]);
            const procid_t source_proc(vertex2proc.at(source));
            const procid_t source_atom(vertex2atomid.at(source));
            const size_t source_idx(source_proc + threadid * rmi.numprocs());
            // add source atom to target
            if(target_proc == rmi.procid()) {
              add_nbr_atom_local(target, source_atom);
            } else {
              proc2buffer.at(target_idx).push_back(args(target, source_atom));
              if(proc2buffer.at(target_idx).size() > BUFFER_SIZE) {
                rmi.remote_call(target_proc,
                                &alist_shuffler_type::add_nbr_atom_local_vec,
                                proc2buffer.at(target_idx));
                proc2buffer.at(target_idx).clear();
              }             
            } // end of add source atom to target
            // add target atom to source
            if(source_proc == rmi.procid()) {
              add_nbr_atom_local(source, target_atom);
            } else {
              proc2buffer.at(source_idx).push_back(args(source, target_atom));
              if(proc2buffer.at(source_idx).size() > BUFFER_SIZE) {
                rmi.remote_call(source_proc,
                                &alist_shuffler_type::add_nbr_atom_local_vec,
                                proc2buffer.at(source_idx));
                proc2buffer.at(source_idx).clear();
              }             
            } // end of add target atom to source
          } // end of for loop
        } // end of parallel for loop 
        
        // Do the extra flush
#pragma omp parallel for num_threads(nthreads)
        for(long i = 0; i < long(proc2buffer.size()); ++i) {
          const procid_t target( i % rmi.numprocs() );
          if(!proc2buffer[i].empty()) {
            rmi.remote_call(target,
                            &alist_shuffler_type::add_nbr_atom_local_vec,
                            proc2buffer[i]);         
          }          
        }
      } // end of compute nbr atoms for all vertices



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
                                &alist_shuffler::add_vertex_local_vec,
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
                            &alist_shuffler::add_vertex_local_vec,
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
      alist_shuffler alist_shuffler(dc);
      dc.full_barrier();
      alist_shuffler.shuffle(path, atom_path);
      dc.full_barrier();
      if(dc.procid() == 0) 
        std::cout << "Finished." << std::endl;      
    } // end of graph partition to atom index






  }; // End of atom shuffler





  












}; // end namespace graphlab














#include <graphlab/macros_undef.hpp>






#endif
