/*
 * Copyright (c) 2014 Daniel McEnnis.
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
 *      https://www.github.com/dmcennis/graphlab
 *
 */

#include <stdlib.h>
#include <graphlab.hpp>

class DjikstraNode {
public:
    unsigned long id;
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

class Gather {
public:
    long id;
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

// The vertex data is its label
typedef PrestigeAnalysisNode vertex_data_type;

typedef GatherMultiTree gather_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, double> graph_type;

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

struct djikstra_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t";
    double value = 0.0;
    for(std::map<long,DjikstraNode>::const_iterator iter = v.data().djikstra_pieces.begin();
        iter != v.data().djikstra_pieces.end();++iter){
        value += iter->second.cost;
    }
    strm << value << std::endl;
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};

size_t num_vertices = 3000;

bool selectVertices(const graph_type::vertex_type& vertex){
    float r = ((float)random()) / ((float)RAND_MAX);
    std::cout << "Random seed is " << r << std::endl;
    if(r < (3000.0 / num_vertices)){
          return true;
    }
    return false;
}



class PrestigeGather{
public:
    std::map<long,double> counts;
    std::map<long,long> edge_count;
	
 void save(graphlab::oarchive& oarc) const {
    oarc << counts << edge_count;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> counts >> edge_count;
  }

  PrestigeGather& operator+=(const PrestigeGather& gather){
    for(std::map<long, double>::const_iterator iter = this->counts.begin();
        iter != this->counts.end(); ++iter ){
                long key = iter->first;
                this->counts[key] += gather.counts.find(key)->second;
                this->edge_count[key] += gather.edge_count.find(key)->second;
	}
    for(std::map<long, double>::const_iterator iter = gather.counts.begin();
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


class PrestigeAlgorithm :
  public graphlab::ivertex_program<graph_type, PrestigeGather>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

    PrestigeGather gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    PrestigeGather g;
    for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
        iter != vertex.data().djikstra_pieces.end(); ++iter){
        long key= iter->first;
        if(edge.target().data().djikstra_pieces[key].id == vertex.id()){
            if(edge.source().data().djikstra_pieces[key].launched == true){
                g.counts[key] = edge.source().data().djikstra_pieces[key].cost + edge.data();
                g.edge_count[key] = 1;
            }
        }
	}
    return g;
    }

    void apply(icontext_type& context, vertex_type& vertex, const PrestigeGather& total) {
        for(std::map<long, DjikstraNode>::const_iterator iter = vertex.data().djikstra_pieces.begin();
                iter != vertex.data().djikstra_pieces.end(); ++iter){
            long key = iter->first;
            if(total.edge_count.find(key)->second==0){
                vertex.data().djikstra_pieces[key].launched = true;
            }
            if((vertex.data().djikstra_pieces[key].launched == true)&&
                    (vertex.data().djikstra_pieces[key].done == false)&&
                    (((long)vertex.data().djikstra_pieces[key].cost)==total.edge_count.find(key)->second)){
                vertex.data().djikstra_pieces[key].done = true;
                vertex.data().djikstra_pieces[key].cost = (double)total.edge_count.find(key)->second;
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
    graphlab::command_line_options clopts("Label Propagation algorithm.");
    std::string graph_dir;
    std::string execution_type = "synchronous";
    clopts.attach_option("graph", graph_dir, "The graph file. Required ");
    clopts.add_positional("graph");
    clopts.attach_option("execution", execution_type, "Execution type (synchronous or asynchronous)");

    std::string saveprefix;
    clopts.attach_option("saveprefix", saveprefix,
                         "If set, will save the resultant pagerank to a "
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

    graphlab::omni_engine<DjikstraAlgorithm> engine(dc, graph, execution_type, clopts);

    num_vertices = graph.num_vertices();
    graphlab::vertex_set start_set = graph.select(selectVertices);
    engine.signal_vset(start_set);
    engine.start();

    const float runtime = engine.elapsed_seconds();
    dc.cout() << "Finished Djikstra engine in " << runtime << " seconds." << std::endl;

    graphlab::omni_engine<ClearBooleans> engine2(dc,graph,execution_type,clopts);
    engine2.signal_all();
    engine2.start();

    graphlab::omni_engine<PrestigeAlgorithm> engine3(dc,graph,execution_type,clopts);
    engine3.signal_all();
    engine3.start();

    if (saveprefix != "") {
      graph.save(saveprefix, djikstra_writer(),
         false,  // do not gzip
         true,   //save vertices
         false); // do not save edges
    }


    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
}


