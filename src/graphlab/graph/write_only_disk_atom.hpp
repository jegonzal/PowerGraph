#ifndef WRITE_ONLY_DISK_ATOM_HPP
#define WRITE_ONLY_DISK_ATOM_HPP

#include <fstream>
#include <string>
#include <graphlab/parallel/parallel_includes.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph_atom.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <graphlab/graph/memory_atom.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

namespace graphlab {
  
  class write_only_disk_atom: public graph_atom {
  private:
    typedef graph<bool,bool>::vertex_id_type    vertex_id_type;
    typedef graph<bool,bool>::vertex_color_type vertex_color_type;
    
    std::string filename;
    uint16_t atomid;
    mutex mut;

    std::ofstream rawofile;
    boost::iostreams::filtering_stream<boost::iostreams::output> f; 

    // some tracking so I can get some basic aggregate reads
    size_t nverts;
    size_t nedges;
    vertex_color_type maxcolor;
    std::map<uint16_t, uint32_t> adjatoms;
    
    inline void open_file(bool clearfile) {
      if (clearfile) {
        rawofile.open(filename.c_str(), std::ios::binary);
      }
      else {
        rawofile.open(filename.c_str(), std::ios::binary | std::ios::app);
      }
      ASSERT_TRUE(rawofile.good());
      f.push(boost::iostreams::zlib_compressor(boost::iostreams::zlib::default_compression, 1024*1024));
      f.push(rawofile);
    }
    
    inline void close_file() {
      f.flush();
      f.pop();
      f.pop();
      rawofile.close();
    }
  public:
    /** constructor. Accesses an atom stored at the filename provided. Clears
     * the existing file if any
    */
    inline write_only_disk_atom(std::string filename, uint16_t atomid, bool clearfile):
                                    filename(filename),atomid(atomid) {
      open_file(false);
      clear_tracking();
    }


    inline void clear_tracking() {
      nverts = 0;
      nedges = 0;
      maxcolor = 0;
      adjatoms.clear();
      
    }
  
    inline ~write_only_disk_atom() { 
      close_file();
    }

    /// Gets the atom ID of this atom
    inline uint16_t atom_id() const {
      return atomid;
    }
  
    inline std::string get_filename() const {
      return filename;
    }
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, it will be overwritten.
     */
    inline void add_vertex(vertex_id_type vid, uint16_t owner) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'a' << vid << owner;
      mut.unlock();
    }
  
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, nothing will be done.
     * Returns true if vertex was added.
     */
    inline bool add_vertex_skip(vertex_id_type vid, uint16_t owner) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'b' << vid << owner;
      mut.unlock();
      return true;
    }
  
  
    inline void add_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string& data) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'c' << vid << owner << data;
      if (owner == atom_id()) nverts++;
      mut.unlock();
    }

  
    inline void add_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) {
      oarchive oarc(f);
      mut.lock();
      if (edata.length() > 0) nedges++;
      oarc << 'f' << src << target << edata;
      mut.unlock();
    }
  
    void add_edge_with_data(vertex_id_type src, uint16_t srcowner,
                           vertex_id_type target, uint16_t targetowner, const std::string &edata) {
      oarchive oarc(f);
      mut.lock();
      if (edata.length() > 0) nedges++;
      if (srcowner != atom_id()) adjatoms[srcowner]++;
      if (targetowner != atom_id()) adjatoms[targetowner]++;
      oarc << 'd' << src << srcowner << target << targetowner <<  edata;
      mut.unlock();
    }

  
    inline void set_vertex(vertex_id_type vid, uint16_t owner) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'g' << vid << owner;
      mut.unlock();
    }
  
    inline void set_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &s) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'h' << vid << owner << s;
      mut.unlock();
    }


    inline void set_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &s) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'j' << src << target << s;
      mut.unlock();
    }
  

  
    inline bool get_vertex(vertex_id_type vid, uint16_t &owner) { ASSERT_TRUE(false); return true;}
  
    inline bool get_vertex_data(vertex_id_type vid, uint16_t &owner, std::string &s) { ASSERT_TRUE(false); return true;}

    inline bool get_edge_data(vertex_id_type src, vertex_id_t target, std::string &s) { ASSERT_TRUE(false); return true;}


    inline std::vector<vertex_id_type> enumerate_vertices() { 
      ASSERT_TRUE(false); 
      return std::vector<vertex_id_type>(); 
    }
  
    /**
     * \brief Returns a list of all the adjacent atoms in the file
     * and the number of ghost vertices in this atom belonging to the
     * adjacent atom
     */
    inline std::map<uint16_t, uint32_t> enumerate_adjacent_atoms() { 
      return adjatoms;
    }
  
    /**
     * \brief Returns the set of incoming vertices of vertex 'vid'
     */
    inline std::vector<vertex_id_type> get_in_vertices(vertex_id_type vid) { 
      ASSERT_TRUE(false); 
      return std::vector<vertex_id_type>(); 
    }
   
   
    /**
     * \brief Returns the set of outgoing vertices of vertex 'vid'
     */
    inline std::vector<vertex_id_type> get_out_vertices(vertex_id_type vid) { 
      ASSERT_TRUE(false); 
      return std::vector<vertex_id_type>(); 
    }


    /**
     * \brief Get the color of the vertex 'vid'.
     * Returns vertex_color_type(-1) if the entry does not exist
     */
    inline vertex_color_type get_color(vertex_id_type vid) { 
      ASSERT_TRUE(false); 
      return vertex_color_type(); 
    }

    /**
     * \brief Sets the color of vertex 'vid'
     */
    inline void set_color(vertex_id_type vid, vertex_color_type color) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'k' << vid << color;
      if (color != vertex_color_type(-1)) {
        maxcolor = std::max(color, maxcolor);
      }
      mut.unlock();
    }
  
    /// Returns the largest color number
    inline vertex_color_type max_color() { return maxcolor; }
    
  
    /**
     * \brief Reads from the auxiliary hash table mapping vid ==> owner.
     * Returns (uint16_t)(-1) if the entry does not exist
     */
    inline uint16_t get_owner(vertex_id_type vid) { ASSERT_TRUE(false); return 0; }

    /**
     * \brief Writes to the auxiliary hash table mapping vid ==> owner.
     */
    inline void set_owner(vertex_id_type vid, uint16_t owner) {
      oarchive oarc(f);
      mut.lock();
      oarc << 'l' << vid << owner;
      mut.unlock();
    }

    /// \brief empties the atom file
    inline void clear() { 
      f.flush();
      rawofile.close();
      rawofile.open(filename.c_str(), std::ios::binary);
      clear_tracking();
    };
  
    /// \brief Ensures the disk storage is up to date
    inline void synchronize() { 
      f.flush();
      rawofile.flush();
    };
    
    /** \brief Return the total number of vertices stored in this atom, 
     * whether or not the this atom actually owns the vertex.
     */
    inline uint64_t num_vertices() const { return nverts; }
  
    /** \brief  Return the total number of edges stored in this atom, whether or 
     * not the this atom actually owns the edge.
     */
    inline uint64_t num_edges() const { return nedges; }
  
    void play_back(graph_atom* atom);
  };

}

#endif

