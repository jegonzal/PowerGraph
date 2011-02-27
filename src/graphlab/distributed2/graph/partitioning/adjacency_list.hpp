#ifndef GRAPHLAB_ADJACENCY_LIST
#define GRAPHLAB_ADJACENCY_LIST


#include <vector>
#include <boost/unordered_map.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/serialization/serialization_includes.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {


  struct adjacency_list {
    //    typedef std::map<vertex_id_t, vertex_id_t> global2local_type;
    typedef boost::unordered_map<vertex_id_t, vertex_id_t> global2local_type;

    static const std::string elist_suffix;
    static const std::string vlist_suffix;
    static const std::string vdata_suffix;
    static const std::string edata_suffix;

    std::vector< vertex_id_t >               local_vertices;
    std::vector< std::vector<vertex_id_t> >  in_nbr_ids;
    global2local_type                        global2local;
    size_t                                   nedges;
   
    //! Load the structure file from an adjacency list a the given
    //! location
    adjacency_list() : nedges(0) {}

    vertex_id_t add_vertex(const vertex_id_t& vid);

    vertex_id_t get_local_vid(const vertex_id_t& gvid) const;

    void add_edge(const vertex_id_t& source, 
                  const vertex_id_t& target,
                  const bool require_target_ownership = false);



    void load(const std::string& fname, const double acceptance_rate = 1);  

    void save(const std::string& base, const size_t& id) const; 

    void operator+=(const adjacency_list& other);

    static std::string 
    make_fname(const std::string& base,
               const size_t& id,
               const std::string& suffix);

    static std::string
    change_suffix(const std::string& fname,
                  const std::string& new_suffix);
    
    static void list_vlist_files(const std::string& pathname, 
                                 std::vector<std::string>& files);
    
    // returns false if local structures are invalid
    void check_local_structures(size_t nverts) const;
  };







}; // end namespace
#include <graphlab/macros_undef.hpp>
#endif





