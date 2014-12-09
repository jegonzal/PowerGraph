/*
 * Copyright (c) 2014 Daniel McEnnis.
 * portions of main Copyright (c) 2009 Carnegie Mellon
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
 *
 */

#include <stdlib.h>
#include <math.h>
#include <graphlab.hpp>

/*
 * Djikstra Graph Node Class
 *
 * This class contains the information about a single graphlab node.
 * id - current best path's previous node id - next node on path to root
 * cost - current cost of the path to route by the current route: Note - this
 *       can become inaccurate in the course of calculations and must be recalculated
 *       by traversing the shortest path tree to get an accurate result.
 * launched - has execution of this node been sheduled
 * done - has execution of this node been completed
 */
class DjikstraNode {
public:
    long id;
    double cost;
    bool launched;
    bool done;

    DjikstraNode(){
        id = 0;
        cost = 1e100;
        launched = false;
        done=false;
    }

  void save(graphlab::oarchive& oarc) const {
    oarc << id << cost << launched << done;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> id >> cost >> launched >> done;
  }
};

/*
 * PrestigeAnalysisNode
 * Graph Node class for running multiple djikstra tree algorithms simultaneously
 * Contains a map of node id's to DjikstraNode instances
 * bookkeeping components
 *
 */
class PrestigeAnalysisNode {
public:
    std::map<long,DjikstraNode> djikstra_pieces;
    double local_value;
    double total;
    long count;
    int edge_count;

    PrestigeAnalysisNode(){
        local_value=0.0;
        total=0.0;
        count=0;
        edge_count=-1;
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << djikstra_pieces << local_value << total << count << edge_count;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> djikstra_pieces >> local_value >> total >> count >> edge_count;
    }
};

/*
 * Gather class for the Djikstra algorithm.
 * id: node id of the incoming edge's other end
 * cost: shortest path cost at the time this node gathers its edges
 * edge_count: a count of gathered edges
 *
 */
class Gather {
public:
    unsigned long id;
    double cost;
    int edge_count;

    Gather(){
        id=0;
        cost=0.0;
        edge_count=1;
    }

    Gather& operator+=(const Gather& other){
        if(other.id < 0){
            return *this;
        }
        if(this->id < 0){
            return *this;
        }
        if (cost <= other.cost){
            this->edge_count++;
            return *this;
        }
        this->edge_count += other.edge_count;
        return *this;
    }


    void save(graphlab::oarchive& oarc) const {
       oarc << id << cost << edge_count;
     }

     void load(graphlab::iarchive& iarc) {
       iarc >> id >> cost >> edge_count;
     }

};

/*
 * GatherMultiTree
 * map of djisktra root id's to their asociated content for that tree
 *
 */
class GatherMultiTree {
public:
    std::map<long,Gather> content;
    int edge_count;

    GatherMultiTree(){
        edge_count=0;
    }

    GatherMultiTree& operator+=(const GatherMultiTree& other){
        return *this;
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << content << edge_count;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> content >> edge_count;
    }
};


typedef PrestigeAnalysisNode vertex_data_type;

typedef GatherMultiTree gather_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, double> graph_type;

/*
 * Loads graphs in the form 'id (id edge_strength)*'
 *
 */
bool line_parser(graph_type& graph, const std::string& filename, const std::string& textline) {
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;
  // first entry in the line is a vertex ID
  strm >> vid;
  PrestigeAnalysisNode node;
  // insert this vertex with its label
  graph.add_vertex(vid, node);
  // while there are elements in the line, continue to read until we fail
  double edge_val=1.0;
  while(1){
    graphlab::vertex_id_type other_vid;
    strm >> other_vid;
    strm >> edge_val;
    if (strm.fail())
      break;
    graph.add_edge(vid, other_vid,edge_val);
  }

  return true;
}

/*
 * Algorithm class whose sole purpose is to reset launched and done booleans
 * for all id's in a PrestigeAnalysisNode
 */
class ClearBooleans :
        public graphlab::ivertex_program<graph_type, gather_type>,
        public graphlab::IS_POD_TYPE {
public:
  edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
        GatherMultiTree g;
        return g;
  }

  void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
      for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
          iter != vertex.data().djikstra_pieces.end(); ++iter){
          long key = iter->first;
          vertex.data().djikstra_pieces[key].launched = false;
          vertex.data().djikstra_pieces[key].done = false;
          vertex.data().djikstra_pieces[key].cost = 0.0;
      }
  }

  edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
          return graphlab::NO_EDGES;
  }

  void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
  }
};

/*
 * Djikstra Algorithm Class
 *
 * Starting from the starting nodes, create an id for this root and signal
 * all neighbors to start the calculations. Set launched when started, done
 * when all edges have been signaled.
 *
 * As a signal is receieved collect edges to determine if the best path has
 * changed.  If it has, update. If the first signal is receieved, marked
 * the node as launched and then mark the node done after signaling neighbors.
 *
 * The process terminates when all nodes active have no neighbors that are not done.
 */
class DjikstraAlgorithm :
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    Gather g;
    GatherMultiTree tree;
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
            iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key=iter->first;
            if((edge.source().data().djikstra_pieces[key].launched == true)&&
                    (edge.source().data().djikstra_pieces[key].done == false)){
                double c = edge.data() + edge.source().data().djikstra_pieces[key].cost;
                g.cost = c;
                g.id = edge.source().data().djikstra_pieces[key].id;
                g.edge_count = 1;
                tree.content[key] = g;
            }else{
                g.id=0;
            }
        }
    return tree;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
            iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key = iter->first;
            if(vertex.data().djikstra_pieces[key].launched == false){
                vertex.data().djikstra_pieces[key].launched = true;
                vertex.data().edge_count = total.edge_count;
                    if(vertex.data().djikstra_pieces[key].cost > total.content.find(key)->second.cost){
                    vertex.data().djikstra_pieces[key].cost = total.content.find(key)->second.cost;
                    vertex.data().djikstra_pieces[key].id = total.content.find(key)->second.id;
                    }else{
                    vertex.data().djikstra_pieces[key].done = true;
                    }
            }else{
                vertex.data().djikstra_pieces[key].done = true;
            }
        }
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
            iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key = iter->first;
            if(vertex.data().djikstra_pieces.find(key)==vertex.data().djikstra_pieces.end()){
                vertex.data().djikstra_pieces[key].launched = true;
                vertex.data().edge_count = total.edge_count;
                vertex.data().djikstra_pieces[key].cost = total.content.find(key)->second.cost;
                vertex.data().djikstra_pieces[key].id = total.content.find(key)->second.id;
            }
        }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      // if vertex data changes, scatter to all edges.
        bool done = true;
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
            iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key = iter->first;
            if(vertex.data().djikstra_pieces.find(key)->second.launched &&
                    !vertex.data().djikstra_pieces.find(key)->second.done){
                done = false;
            }
        }
        if(!done){
            return graphlab::OUT_EDGES;
        }else{
            return graphlab::NO_EDGES;
        }
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
        iter != vertex.data().djikstra_pieces.end(); ++iter){
        long key = iter->first;
        if((vertex.data().djikstra_pieces.find(key)->second.done == false) &&
                (vertex.data().djikstra_pieces.find(key)->second.launched == true)){
                context.signal(edge.target());
            }
    }
  }
};

size_t num_vertices = 3000;
size_t desired_vertices_count = 3000;
size_t selected_vertices_count = 0;

/*
 * For every node, print the previous node in its spanning tree for all spanning trees this node is in.
 *
 */
struct betweeness_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id();
    double betweeness = 0.0;
    for(std::map<long, DjikstraNode>::const_iterator iter = v.data().djikstra_pieces.begin();
        iter != v.data().djikstra_pieces.end(); ++iter){
        betweeness += iter->second.cost;
    }
    betweeness /= selected_vertices_count;
    strm << "\t" << betweeness << std::endl;
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};

/*
 * Select ~3000 root nodes or an exact count which gives up around +/-3% accuracy
 * in prestige measures. It is a constant memory random selector.
 */
bool selectVertices(const graph_type::vertex_type& vertex){
    unsigned int r = random();
    std::cout << "Random seed is " << r << std::endl;
    if(r < (desired_vertices_count * RAND_MAX / num_vertices)){
          selected_vertices_count++;
          return true;
    }
    return false;
}


/*
 * Gather object that keeps track of betweeness counts for each spanning tree.
 *
 */
class BetweenessGather{
public:
    std::map<long,long> counts;
    std::map<long,long> edge_count;
	
 void save(graphlab::oarchive& oarc) const {
    oarc << counts << edge_count;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> counts >> edge_count;
  }

  BetweenessGather& operator+=(const BetweenessGather& gather){
    for(std::map<long, long>::const_iterator iter = this->counts.begin();
        iter != this->counts.end(); ++iter ){
                long key = iter->first;
                this->counts[key] += gather.counts.find(key)->second;
                this->edge_count[key] += gather.edge_count.find(key)->second;
	}
    for(std::map<long, long>::const_iterator iter = gather.counts.begin();
            iter != gather.counts.end(); ++iter){
                long key = iter->first;
                if(this->counts.find(key) != this->counts.end()){
                        this->counts[key] = gather.counts.find(key)->second;
                        this->edge_count[key] = gather.edge_count.find(key)->second;
                }
    }
    return *this;
  }


};

/*
 * Walk backwards from leaf nodes (those that have no nodes pointing to them in
 * the gather step). Each signals the node referenced in its internal spanning tree
 * record. This is performed simultaneously for each spanning tree in the set.
 *
 * The betweeness score is cached in the cost field.
 *
 */
class BetweenessAlgorithm :
  public graphlab::ivertex_program<graph_type, BetweenessGather>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

    BetweenessGather gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
	BetweenessGather g;
    for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
        iter != vertex.data().djikstra_pieces.end(); ++iter){
        long key= iter->first;
        if(edge.target().data().djikstra_pieces[key].id == vertex.id()){
            if(edge.source().data().djikstra_pieces[key].launched == true){
                g.counts[key] = edge.source().data().djikstra_pieces[key].cost;
                g.edge_count[key] = 1;
            }
        }
	}
    return g;
    }

    void apply(icontext_type& context, vertex_type& vertex, const BetweenessGather& total) {
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
                iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key = iter->first;
            if(total.edge_count.find(key)->second==0){
                vertex.data().djikstra_pieces[key].launched = true;
                vertex.data().djikstra_pieces[key].cost = 0.0;
            }
            if((vertex.data().djikstra_pieces[key].launched == true)&&
                    (vertex.data().djikstra_pieces[key].done == false)&&
                    (((long)vertex.data().djikstra_pieces[key].cost)==total.edge_count.find(key)->second)){
                vertex.data().djikstra_pieces[key].done = true;
                vertex.data().djikstra_pieces[key].cost = fmax(1.0,(double)total.edge_count.find(key)->second);
            }
        }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      // if vertex data changes, scatter to all edges.
     	bool done = true;
    for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
        iter != vertex.data().djikstra_pieces.end(); ++iter){
        long key = iter->first;
        if(vertex.data().djikstra_pieces.find(key)->second.launched && !vertex.data().djikstra_pieces.find(key)->second.done){
        		done = false;
      		}
        }
        if(!done){
                return graphlab::OUT_EDGES;
        }else{
                return graphlab::NO_EDGES;
        }
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
            iter != vertex.data().djikstra_pieces.end(); ++iter){
                long key = iter->first;
                if((vertex.data().djikstra_pieces.find(key)->second.done == false) &&
                        (vertex.data().djikstra_pieces.find(key)->second.launched == true)){
                        context.signal(edge.target());
                }
        }
    }
};

int main (int argc, char** argv){
    // Initialize control plain using mpi
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
    global_logger().set_log_level(LOG_INFO);

    // Parse command line options -----------------------------------------------
    graphlab::command_line_options clopts("Betweeness Algorithm.");
    std::string graph_dir;
    clopts.attach_option("graph", graph_dir, "The graph file. Required ");
    clopts.add_positional("graph");
    clopts.attach_option("samplesize", desired_vertices_count, "(Sample Size) Number of spanning trees to use");

    std::string saveprefix;
    clopts.attach_option("saveprefix", saveprefix,
                         "If set, will save the resultant betweness score to a "
                         "sequence of files with prefix saveprefix");

    if(!clopts.parse(argc, argv)) {
      dc.cout() << "Error in parsing command line arguments." << std::endl;
      return EXIT_FAILURE;
    }
    if (graph_dir == "") {
      dc.cout() << "Graph not specified. Cannot continue";
      return EXIT_FAILURE;
    }

    // Build the graph ----------------------------------------------------------
    graph_type graph(dc);
    dc.cout() << "Loading graph using line parser" << std::endl;
    graph.load(graph_dir, line_parser);

    dc.cout() << "#vertices: " << graph.num_vertices() << " #edges:" << graph.num_edges() << std::endl;

    graphlab::omni_engine<DjikstraAlgorithm> engine(dc, graph, "asynchronous", clopts);

    num_vertices = graph.num_vertices();
    graphlab::vertex_set start_set = graph.select(selectVertices);
    engine.signal_vset(start_set);
    engine.start();

    const float runtime = engine.elapsed_seconds();
    dc.cout() << "Finished Djikstra engine in " << runtime << " seconds." << std::endl;

    graphlab::omni_engine<ClearBooleans> engine2(dc,graph,"asynchronous",clopts);
    engine2.signal_all();
    engine2.start();
    const float runtime2 = engine.elapsed_seconds();
    dc.cout() << "Finished resetting graph engine in " << runtime2 << " seconds." << std::endl;

    graphlab::omni_engine<BetweenessAlgorithm> engine3(dc,graph,"asynchronous",clopts);
    engine3.signal_all();
    engine3.start();
    const float runtime3 = engine.elapsed_seconds();
    dc.cout() << "Finished Betweeness engine in " << runtime3 << " seconds." << std::endl;

    if (saveprefix != "") {
      graph.save(saveprefix, betweeness_writer(),
         false,  // do not gzip
         true,   //save vertices
         false); // do not save edges
    }


    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
}


