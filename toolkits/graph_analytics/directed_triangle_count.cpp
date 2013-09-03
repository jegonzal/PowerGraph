/*  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#include <boost/unordered_set.hpp>
#include <graphlab.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/util/cuckoo_set_pow2.hpp>
#include <graphlab/macros_def.hpp>
/**
 This implements the exact counting procedure described in 

 Efficient Algorithms for Large-Scale Local Triangle Counting
 Luca Becchetti, Paolo Boldi, Carlos Castillo, Aristides Gioni  

  */
   

// Radix sort implementation from https://github.com/gorset/radix
// Thanks to Erik Gorset
//
/*
Copyright 2011 Erik Gorset. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Erik Gorset ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Erik Gorset OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Erik Gorset.
*/
void radix_sort(graphlab::vertex_id_type *array, int offset, int end, int shift) {
    int x, y;
    graphlab::vertex_id_type value, temp;
    int last[256] = { 0 }, pointer[256];

    for (x=offset; x<end; ++x) {
        ++last[(array[x] >> shift) & 0xFF];
    }

    last[0] += offset;
    pointer[0] = offset;
    for (x=1; x<256; ++x) {
        pointer[x] = last[x-1];
        last[x] += last[x-1];
    }

    for (x=0; x<256; ++x) {
        while (pointer[x] != last[x]) {
            value = array[pointer[x]];
            y = (value >> shift) & 0xFF;
            while (x != y) {
                temp = array[pointer[y]];
                array[pointer[y]++] = value;
                value = temp;
                y = (value >> shift) & 0xFF;
            }
            array[pointer[x]++] = value;
        }
    }

    if (shift > 0) {
        shift -= 8;
        for (x=0; x<256; ++x) {
            temp = x > 0 ? pointer[x] - pointer[x-1] : pointer[0] - offset;
            if (temp > 64) {
                radix_sort(array, pointer[x] - temp, pointer[x], shift);
            } else if (temp > 1) {
                std::sort(array + (pointer[x] - temp), array + pointer[x]);
                //insertion_sort(array, pointer[x] - temp, pointer[x]);
            }
        }
    }
}

size_t HASH_THRESHOLD = 64;

// We on each vertex, either a vector of sorted VIDs
// or a hash set (cuckoo hash) of VIDs.
// If the number of elements is greater than HASH_THRESHOLD,
// the hash set is used. Otherwise the vector is used.
struct vid_vector{
  std::vector<graphlab::vertex_id_type> vid_vec;
  graphlab::cuckoo_set_pow2<graphlab::vertex_id_type, 3> *cset;
  vid_vector(): cset(NULL) { }
  vid_vector(const vid_vector& v):cset(NULL) {
    (*this) = v;
  }

  vid_vector& operator=(const vid_vector& v) {
    if (this == &v) return *this;
    vid_vec = v.vid_vec;
    if (v.cset != NULL) {
      // allocate the cuckoo set if the other side is using a cuckoo set
      // or clear if I alrady have one
      if (cset == NULL) {
        cset = new graphlab::cuckoo_set_pow2<graphlab::vertex_id_type, 3>(-1, 0, 2 * v.cset->size());
      }
      else {
        cset->clear();
      }
      (*cset) = *(v.cset);
    }
    else {
      // if the other side is not using a cuckoo set, lets not use a cuckoo set
      // either
      if (cset != NULL) {
        delete cset;
        cset = NULL;
      }
    }
    return *this;
  }

  ~vid_vector() {
    if (cset != NULL) delete cset;
  }

  // assigns a vector of vertex IDs to this storage.
  // this function will clear the contents of the vid_vector
  // and reconstruct it.
  // If the assigned values has length >= HASH_THRESHOLD,
  // we will allocate a cuckoo set to store it. Otherwise,
  // we just store a sorted vector
  void assign(const std::vector<graphlab::vertex_id_type>& vec) {
    clear();
    if (vec.size() >= HASH_THRESHOLD) {
        // move to cset
        cset = new graphlab::cuckoo_set_pow2<graphlab::vertex_id_type, 3>(-1, 0, 2 * vec.size());
        foreach (graphlab::vertex_id_type v, vec) {
          cset->insert(v);
        }
    }
    else {
      vid_vec = vec;
      if (vid_vec.size() > 64) {
        radix_sort(&(vid_vec[0]), 0, vid_vec.size(), 24);
      }
      else {
        std::sort(vid_vec.begin(), vid_vec.end());
      }
      std::vector<graphlab::vertex_id_type>::iterator new_end = std::unique(vid_vec.begin(),
                                               vid_vec.end());
      vid_vec.erase(new_end, vid_vec.end());
    }
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << (cset != NULL);
    if (cset == NULL) oarc << vid_vec;
    else oarc << (*cset);
  }


  void clear() {
    vid_vec.clear();
    if (cset != NULL) {
      delete cset;
      cset = NULL;
    }
  }

  size_t size() const {
    return cset == NULL ? vid_vec.size() : cset->size();
  }

  void load(graphlab::iarchive& iarc) {
    clear();
    bool hascset;
    iarc >> hascset;
    if (!hascset) iarc >> vid_vec;
    else {
      cset = new graphlab::cuckoo_set_pow2<graphlab::vertex_id_type, 3>(-1, 0, 2);
      iarc >> (*cset);
    }
  }
};

/*
  A simple counting iterator which can be used as an insert iterator.
  but only counts the number of elements inserted. Useful for
  use with counting the size of an intersection using std::set_intersection
*/
template <typename T>
struct counting_inserter {
  size_t* i;
  counting_inserter(size_t* i):i(i) { }
  counting_inserter& operator++() {
    ++(*i);
    return *this;
  }
  void operator++(int) {
    ++(*i);
  }

  struct empty_val {
    empty_val operator=(const T&) { return empty_val(); }
  };

  empty_val operator*() {
    return empty_val();
  }

  typedef empty_val reference;
};


/*
 * Computes the size of the intersection of two vid_vector's
 */
static uint32_t count_set_intersect(
             const vid_vector& smaller_set,
             const vid_vector& larger_set) {
  if (smaller_set.size() > larger_set.size()) {
    return count_set_intersect(larger_set, smaller_set);
  }
  if (smaller_set.cset == NULL && larger_set.cset == NULL) {
    size_t i = 0;
    counting_inserter<graphlab::vertex_id_type> iter(&i);
    std::set_intersection(smaller_set.vid_vec.begin(), smaller_set.vid_vec.end(),
                          larger_set.vid_vec.begin(), larger_set.vid_vec.end(),
                          iter);
    return i;
  }
  else if (smaller_set.cset == NULL && larger_set.cset != NULL) {
    size_t i = 0;
    foreach(graphlab::vertex_id_type vid, smaller_set.vid_vec) {
      i += larger_set.cset->count(vid);
    }
    return i;
  }
  else if (smaller_set.cset != NULL && larger_set.cset == NULL) {
    size_t i = 0;
    foreach(graphlab::vertex_id_type vid, larger_set.vid_vec) {
      i += smaller_set.cset->count(vid);
    }
    return i;
  }
  else {
    size_t i = 0;
    foreach(graphlab::vertex_id_type vid, *(smaller_set.cset)) {
      i += larger_set.cset->count(vid);
    }
    return i;

  }
}


// This structure is used to hold the final triangle counts 
// on each vertex
struct triangle_count: public graphlab::IS_POD_TYPE {
  triangle_count(): out_triangles(0), in_triangles(0), 
                      through_triangles(0), cycle_triangles(0) { }
  // A is the example below
  /*
    A ---> B
    |   / 
    |  /
    v /
    C
           diagonal edge direction does not matter
  */
  uint32_t out_triangles;

/*
    A<--- B
    ^   / 
    |  /
    | /
    C
          diagonal edge direction does not matter
  */
  uint32_t in_triangles;

  /*
    A---> B
    ^   ^ 
    |  /
    | /
    C
  */
  uint32_t through_triangles;

  /*
    A---> B
    ^   / 
    |  /
    | v
    C
  */
  uint32_t cycle_triangles;

  triangle_count& operator+=(const triangle_count& other) {
    out_triangles += other.out_triangles;
    in_triangles += other.in_triangles;
    through_triangles += other.through_triangles;
    cycle_triangles += other.cycle_triangles;
    return *this;
  }
};

/*
 * Each vertex maintains a list of all its neighbors.
 * and a final count for the number of triangles it is involved in
 */
struct vertex_data_type {
  vertex_data_type(){ }
  // A list of all its neighbors
  vid_vector in_vid_set;
  vid_vector out_vid_set;
  triangle_count count;
 // The number of triangles this vertex is involved it.
  // only used if "per vertex counting" is used
  void save(graphlab::oarchive &oarc) const {
    oarc << in_vid_set << out_vid_set << count;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> in_vid_set >> out_vid_set >> count;
  }
};




// This structure is used to hold the final triangle counts 
// on each edge
struct edge_triangle_count: public graphlab::IS_POD_TYPE {
  edge_triangle_count(): s_s(0), st_st(0), 
                      st_s(0) { }
  // using notation from the paper
  // s_s is the intersection between outgoing of source and outgoing of target
  // st_st is intersection between incoming of source and incoming of target
  // st_s is intersectino between incoming of source and outgoing of target
  uint32_t s_s, st_st, st_s;

  edge_triangle_count & operator+=(const edge_triangle_count& other) {
    s_s += other.s_s;
    st_st += other.st_st;
    st_s += other.st_s;
    return *this;
  }
};














/*
 * Each edge is simply a counter of triangles
 */
typedef edge_triangle_count edge_data_type;


/*
 * This is the gathering type which accumulates an array of
 * all neighboring vertices.
 * It is a simple wrapper around a vector with
 * an operator+= which simply performs a  +=
 */
struct set_union_gather {
  graphlab::vertex_id_type v;
  std::vector<graphlab::vertex_id_type> vid_vec;

  set_union_gather():v(-1) {
  }

  size_t size() const {
    if (v == (graphlab::vertex_id_type)-1) return vid_vec.size();
    else return 1;
  }
  /*
   * Combining with another collection of vertices.
   * Union it into the current set.
   */
  set_union_gather& operator+=(const set_union_gather& other) {
    if (size() == 0) {
      (*this) = other;
      return (*this);
    }
    else if (other.size() == 0) {
      return *this;
    }

    if (vid_vec.size() == 0) {
      vid_vec.push_back(v);
      v = (graphlab::vertex_id_type)(-1);
    }
    if (other.vid_vec.size() > 0) {
      size_t ct = vid_vec.size();
      vid_vec.resize(vid_vec.size() + other.vid_vec.size());
      for (size_t i = 0; i < other.vid_vec.size(); ++i) {
        vid_vec[ct + i] = other.vid_vec[i];
      }
    }
    else if (other.v != (graphlab::vertex_id_type)-1) {
      vid_vec.push_back(other.v);
    }
    return *this;
  }
  
  // serialize
  void save(graphlab::oarchive& oarc) const {
    oarc << bool(vid_vec.size() == 0);
    if (vid_vec.size() == 0) oarc << v;
    else oarc << vid_vec;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
    bool novvec;
    v = (graphlab::vertex_id_type)(-1);
    vid_vec.clear();
    iarc >> novvec;
    if (novvec) iarc >> v;
    else iarc >> vid_vec;
  }
};


struct set_union_pair{
  set_union_gather in_set;
  set_union_gather out_set;

  set_union_pair& operator+=(const set_union_pair& other) {
    in_set += other.in_set;
    out_set += other.out_set;
    return (*this);
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << in_set << out_set;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
    iarc >> in_set >> out_set;
  }
 
};

/*
 * Define the type of the graph
 */
typedef graphlab::distributed_graph<vertex_data_type,
                                    edge_data_type> graph_type;



/*
 * This class implements the triangle counting algorithm as described in
 * the header. On gather, we accumulate a set of all adjacent vertices.
 * If per_vertex output is not necessary, we can use the optimization
 * where each vertex only accumulates neighbors with greater vertex IDs.
 */
class triangle_count_program :
      public graphlab::ivertex_program<graph_type,
                                      set_union_pair>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE  {
public:
  bool do_not_scatter;

  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } 

  /*
   * For each edge, figure out the ID of the "other" vertex
   * and accumulate a set of the neighborhood vertex IDs.
   */
  gather_type gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    set_union_pair gather;
    // check the edge direction
    if (edge.source().id() == vertex.id()) {
      // this is an out_edge
      graphlab::vertex_id_type otherid = edge.target().id();
      gather.out_set.v = otherid;
    }
    else {
      // this is an in_edge
      graphlab::vertex_id_type otherid = edge.source().id();
      gather.in_set.v = otherid;
    }
    return gather;
  }

  /*
   * the gather result now contains the vertex IDs in the neighborhood.
   * store it on the vertex. 
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& neighborhood) {
   do_not_scatter = false;
   if (neighborhood.in_set.vid_vec.size() == 0) {
     // neighborhood set may be empty or has only 1 element
     vertex.data().in_vid_set.clear();
     if (neighborhood.in_set.v != (graphlab::vertex_id_type(-1))) {
       vertex.data().in_vid_set.vid_vec.push_back(neighborhood.in_set.v);
     }
   }
   else {
     vertex.data().in_vid_set.assign(neighborhood.in_set.vid_vec);
   }


   if (neighborhood.out_set.vid_vec.size() == 0) {
     // neighborhood set may be empty or has only 1 element
     vertex.data().out_vid_set.clear();
     if (neighborhood.out_set.v != (graphlab::vertex_id_type(-1))) {
       vertex.data().out_vid_set.vid_vec.push_back(neighborhood.out_set.v);
     }
   }
   else {
     vertex.data().out_vid_set.assign(neighborhood.out_set.vid_vec);
   }

   do_not_scatter = vertex.data().in_vid_set.size() == 0 && 
                    vertex.data().out_vid_set.size() == 0 ;
  } // end of apply

  /*
   * Scatter over all edges to compute the intersection.
   * I only need to touch each edge once, so if I scatter just on the
   * out edges, that is sufficient.
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    if (do_not_scatter) return graphlab::NO_EDGES;
    else return graphlab::OUT_EDGES;
  }


  /*
   * For each edge, count the intersection of the neighborhood of the
   * adjacent vertices. This is the number of triangles this edge is involved
   * in.
   */
  void scatter(icontext_type& context,
              const vertex_type& vertex,
              edge_type& edge) const {

    //vertex_type othervtx = edge.target();
    // ok. the real work happens here.
    
    const vertex_data_type& srclist = edge.source().data();
    const vertex_data_type& targetlist = edge.target().data();

    edge.data().s_s = count_set_intersect(srclist.out_vid_set,
                                    targetlist.out_vid_set);
    edge.data().st_st += count_set_intersect(srclist.in_vid_set,
                                    targetlist.in_vid_set);
    edge.data().st_s += count_set_intersect(srclist.in_vid_set,
                                    targetlist.out_vid_set);

  }
};

/*
 * This class is used in a second engine call if per vertex counts are needed.
 * The number of triangles a vertex is involved in can be computed easily
 * by summing over the number of triangles each adjacent edge is involved in
 * and dividing by 2. 
 */
class get_per_vertex_count :
      public graphlab::ivertex_program<graph_type, triangle_count>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE  {
public:
  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }
  // We gather the number of triangles each edge is involved in
  triangle_count gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    triangle_count ret;
    if (edge.source().id() == vertex.id()) {
      ret.out_triangles += edge.data().s_s;
      ret.through_triangles += edge.data().st_st;
      ret.cycle_triangles += edge.data().st_s;
    }
    else {
      ret.in_triangles += edge.data().st_st;
    }
    return ret;
  }

  /* the gather result is the total sum of the number of triangles
   * each adjacent edge is involved in . Dividing by 2 gives the
   * desired result.
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& tc) {
    vertex.data().in_vid_set.clear();
    vertex.data().out_vid_set.clear();
    vertex.data().count = tc;
  }

  // No scatter
  edge_dir_type scatter_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }


};

typedef graphlab::synchronous_engine<triangle_count_program> engine_type;

/* Used to sum over all the vertices in the graph in a
 * map_reduce_vertices call
 * to get the total number of triangles
 */
triangle_count get_vertex_counts(const graph_type::vertex_type& v) {
  return v.data().count;
}

/*
 * A saver which saves a file where each line is a vid / # triangles pair
 */
struct save_triangle_count{
  std::string save_vertex(graph_type::vertex_type v) { 
    triangle_count tc = v.data().count;
    double n_followed = v.num_out_edges();
    double n_following = v.num_in_edges();

    return graphlab::tostr(v.id()) + "\t" +
           graphlab::tostr(tc.in_triangles) + "\t" +
           graphlab::tostr(tc.out_triangles) + "\t" +
           graphlab::tostr(tc.through_triangles) + "\t" +
           graphlab::tostr(tc.cycle_triangles) + "\t" +
           graphlab::tostr(n_followed) + "\t" + 
           graphlab::tostr(n_following) + "\n";
  }
  std::string save_edge(graph_type::edge_type e) {
    return "";
  }
};


int main(int argc, char** argv) {
  std::cout << "This program counts the exact number of triangles in the "
            "provided graph.\n\n";

  graphlab::command_line_options clopts("Exact Triangle Counting. "
    "Given a graph, this program computes the total number of triangles "
    "in the graph. An option (per_vertex) is also provided which "
    "computes for each vertex, the number of triangles it is involved in."
    "The algorithm assumes that each undirected edge appears exactly once "
    "in the graph input. If edges may appear more than once, this procedure "
    "will over count.");
  std::string prefix, format;
  std::string per_vertex;
  bool PER_VERTEX_COUNT = false;
  clopts.attach_option("graph", prefix,
                       "Graph input. reads all graphs matching prefix*");
  clopts.attach_option("format", format,
                       "The graph format");
  clopts.attach_option("ht", HASH_THRESHOLD,
                       "Above this size, hash tables are used");
  clopts.attach_option("per_vertex", per_vertex,
                       "If not empty, will count the number of "
                       "triangles each vertex belongs to and "
                       "save to file with prefix \"[per_vertex]\". "
                       "The algorithm used is slightly different "
                       "and thus will be a little slower");
  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (prefix == "") {
    std::cout << "--graph is not optional\n";
    clopts.print_description();
    return EXIT_FAILURE;
  }
  else if (format == "") {
    std::cout << "--format is not optional\n";
    clopts.print_description();
    return EXIT_FAILURE;
  }

  if (per_vertex != "") PER_VERTEX_COUNT = true;
  // Initialize control plane using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  graphlab::launch_metric_server();
  // load graph
  graph_type graph(dc, clopts);
  graph.load_format(prefix, format);
  graph.finalize();
  dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
            << "Number of edges:    " << graph.num_edges() << std::endl;

  graphlab::timer ti;
  
  // create engine to count the number of triangles
  dc.cout() << "Counting Triangles..." << std::endl;
  engine_type engine(dc, graph, clopts);
  engine.signal_all();
  engine.start();

  dc.cout() << "Counted in " << ti.current_time() << " seconds" << std::endl;
  dc.cout() << "Collecting results ... " << std::endl;
  graphlab::synchronous_engine<get_per_vertex_count> engine2(dc, graph, clopts);
  engine2.signal_all();
  engine2.start();
 
  if (PER_VERTEX_COUNT == false) {
    triangle_count count = graph.map_reduce_vertices<triangle_count>(get_vertex_counts);
    dc.cout() << count.in_triangles << " In triangles\n";
    dc.cout() << count.out_triangles << " Out triangles\n";
    dc.cout() << count.through_triangles << " Through triangles\n";
    dc.cout() << count.cycle_triangles << " Cycle triangles\n";
  }
  else {
   dc.cout() << "Saving Results...\n";
   dc.cout() << "Format is \n";
   dc.cout() << "   [vid]  [in triangles]  [out triangles]   [through triangles]  [cycle_triangles]  [#out edges] [#in edges]" << std::endl;
   graph.save(per_vertex,
            save_triangle_count(),
            false, /* no compression */
            true, /* save vertex */
            false, /* do not save edge */
            clopts.get_ncpus()); /* one file per machine */

  }
  
  graphlab::stop_metric_server();

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main

