#ifndef GRAPHLAB_ADJACENCY_LIST
#define GRAPHLAB_ADJACENCY_LIST


#include <vector>
#include <graphlab/graph/graph.hpp>
#include <graphlab/serialization/serialization_includes.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {


  struct adjacency_list {
    typedef std::map<vertex_id_t, vertex_id_t> global2local_type;

    static const std::string alist_suffix;
    static const std::string vdata_suffix;
    static const std::string edata_suffix;

    std::vector< vertex_id_t >               local_vertices;
    std::vector< std::vector<vertex_id_t> >  in_neighbor_ids;
    global2local_type                        global2local;
   
    //! Load the structure file from an adjacency list a the given
    //! location
    adjacency_list(){}

    vertex_id_t add_vertex(const vertex_id_t& vid);

    void add_edge(const vertex_id_t& source, 
                  const vertex_id_t& target);



    void load(const std::string& fname);  

    void save(const std::string& fname) const; 

    void operator+=(const adjacency_list& other);

    static std::string 
    make_fname(const std::string& base,
               const size_t& id,
               const std::string& suffix);

    static std::string 
    make_alist_fname(const std::string& base,
                     const size_t& id);

    static std::string 
    make_vdata_fname(const std::string& base,
                     const size_t& id);

    static std::string 
    make_edata_fname(const std::string& base,
                     const size_t& id);

    static std::string
    change_suffix(const std::string& fname,
                  const std::string& new_suffix);
    
    static void list_adjacency_files(const std::string& pathname, 
                                     std::vector<std::string>& files);
    
  };







}; // end namespace
#include <graphlab/macros_undef.hpp>
#endif





