#ifndef GRAPH_ATOM_HPP
#define GRAPH_ATOM_HPP

#include <sstream>
#include <map>
#include <boost/unordered_map.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>

namespace graphlab {
  
  class graph_atom {
  private:
    typedef graph<bool,bool>::vertex_id_type    vertex_id_type;
    typedef graph<bool,bool>::vertex_color_type vertex_color_type;
  public:
    graph_atom() { };
    virtual ~graph_atom() { };
  
    /// Gets the atom ID of this atom
    virtual uint16_t atom_id() const = 0;
  
    virtual std::string get_filename() const = 0;
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, it will be overwritten.
     */
    virtual void add_vertex(vertex_id_type vid, uint16_t owner) = 0;
  
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, nothing will be done.
     * Returns true if vertex was added.
     */
    virtual bool add_vertex_skip(vertex_id_type vid, uint16_t owner) = 0;
  
  
    virtual void add_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string& data) = 0;
    /**
     * \brief Inserts vertex 'vid' into the file. If the vertex already exists,
     * it will be overwritten.
     */
    template <typename T>
    void add_vertex(vertex_id_type vid, uint16_t owner, const T &vdata) {
      add_vertex_with_data(vid, owner, serialize_to_string(vdata));
    }
  
    /**
     * \brief Inserts edge src->target into the file with data. 
     * If the edge already exists, the data will be overwritten
     */
    virtual void add_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &edata) = 0;
    
    virtual void add_edge_with_data(vertex_id_type src, uint16_t srcowner,
                                    vertex_id_type target, uint16_t targetowner, const std::string &edata) = 0;
    
    
    inline void add_edge(vertex_id_type src, vertex_id_type target) {
      add_edge_with_data(src, target, "");
    }
    
    inline void add_edge(vertex_id_type src, uint16_t srcowner,
                         vertex_id_type target, uint16_t targetowner) {
      add_edge_with_data(src, srcowner, target, targetowner, "");
    }

    /**
     * \brief Inserts edge src->target into the file. If the edge already exists,
     * it will be overwritten.
     */
    template <typename T>
    void add_edge(vertex_id_type src, vertex_id_type target, const T &edata) {
      add_edge_with_data(src, target, serialize_to_string(edata));
    }

    template <typename T>
    void add_edge(vertex_id_type src, uint16_t srcowner,
                 vertex_id_type target, uint16_t targetowner, const T &edata) {
      add_edge_with_data(src, srcowner, target, targetowner, serialize_to_string(edata));
    }
  
  
    /**
     * \brief Modifies an existing vertex in the file where no data is assigned to the 
     * vertex. User must ensure that the file
     * already contains this vertex. If user is unsure, add_vertex should be used.
     */
    virtual void set_vertex(vertex_id_type vid, uint16_t owner) = 0;
  
    virtual void set_vertex_with_data(vertex_id_type vid, uint16_t owner, const std::string &s) = 0;

    /**
     * \brief Modifies an existing vertex in the file. User must ensure that the file
     * already contains this vertex. If user is unsure, add_vertex should be used.
     */
    template <typename T>
    void set_vertex(vertex_id_type vid, uint16_t owner, const T &vdata) {
      set_vertex_with_data(vid, owner, serialize_to_string(vdata));
    }
  
    virtual void set_edge_with_data(vertex_id_type src, vertex_id_type target, const std::string &s) = 0;
  
    template <typename T>
    void set_edge(vertex_id_type src, vertex_id_type target, const T &edata) {
      set_edge_with_data(src, target, serialize_to_string(edata));
    }
  
  
    virtual bool get_vertex(vertex_id_type vid, uint16_t &owner) = 0;
  
    virtual bool get_vertex_data(vertex_id_type vid, uint16_t &owner, std::string &s) = 0;

    /**
     * \brief Reads a vertex from the file returning results in 'owner' and 'vdata'.
     * Returns true if vertex exists and false otherwise.
     * If there is no vertex data stored, vdata will not be modified.
     */
    template <typename T>
    bool get_vertex(vertex_id_type vid, uint16_t &owner, T &vdata) {
      std::string s;
      bool ret = get_vertex_data(vid, owner, s);
      deserialize_from_string(s, vdata);
      return ret;
    }
  
    virtual bool get_edge_data(vertex_id_type src, vertex_id_t target, std::string &s) = 0;

    /**
     * \brief Reads a edge from the file returning results in 'owner' and 'vdata'.
     * Returns true if edge exists and false otherwise.
     * If there is no edge data stored, edata will not be modified.
     */
    template <typename T>
    bool get_edge(vertex_id_type src, vertex_id_type target, T &edata) {
      std::string s;
      bool ret = get_edge_data(src, target, s);
      deserialize_from_string(s, edata);
      return ret;
    }


    virtual std::vector<vertex_id_type> enumerate_vertices() = 0;
  
    /**
     * \brief Returns a list of all the adjacent atoms in the file
     * and the number of ghost vertices in this atom belonging to the
     * adjacent atom
     */
    virtual std::map<uint16_t, uint32_t> enumerate_adjacent_atoms() = 0;
  
    /**
     * \brief Returns the set of incoming vertices of vertex 'vid'
     */
    virtual std::vector<vertex_id_type> get_in_vertices(vertex_id_type vid) = 0;
   
   
    /**
     * \brief Returns the set of outgoing vertices of vertex 'vid'
     */
    virtual std::vector<vertex_id_type> get_out_vertices(vertex_id_type vid) = 0;


    /**
     * \brief Get the color of the vertex 'vid'.
     * Returns vertex_color_type(-1) if the entry does not exist
     */
    virtual vertex_color_type get_color(vertex_id_type vid) = 0;

    /**
     * \brief Sets the color of vertex 'vid'
     */
    virtual void set_color(vertex_id_type vid, vertex_color_type color) = 0;
  
    /// Returns the largest color number
    virtual vertex_color_type max_color() = 0;
  
    /**
     * \brief Reads from the auxiliary hash table mapping vid ==> owner.
     * Returns (uint16_t)(-1) if the entry does not exist
     */
    virtual uint16_t get_owner(vertex_id_type vid) = 0;

    /**
     * \brief Writes to the auxiliary hash table mapping vid ==> owner.
     */
    virtual void set_owner(vertex_id_type vid, uint16_t owner) = 0;

    /// \brief empties the atom file
    virtual void clear() = 0;
  
    /// \brief Ensures the disk storage is up to date
    virtual void synchronize() = 0;
    
    /** \brief Return the total number of vertices stored in this atom, 
     * whether or not the this atom actually owns the vertex.
     */
    virtual uint64_t num_vertices() const = 0;
  
    /** \brief  Return the total number of edges stored in this atom, whether or 
     * not the this atom actually owns the edge.
     */
    virtual uint64_t num_edges() const = 0;
  
  };

}

#endif
