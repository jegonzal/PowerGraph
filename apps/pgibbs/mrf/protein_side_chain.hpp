#ifndef PROTEIN_SIDE_CHAIN_HPP
#define PROTEIN_SIDE_CHAIN_HPP





#include <iostream>
#include <fstream>

#include <vector>
#include <map>

#include <graphlab.hpp>

#include "data_structures.hpp"


/** An object representing the protein data */
struct protein_data {
  std::string name;
  std::vector<unary_factor> potentials;
  std::vector< std::pair<size_t, size_t> > edges;
  std::vector<binary_factor> edge_factors;  
};

std::ostream& operator<<(std::ostream& out, const protein_data& protein) {
  out << "Protein data file: " << protein.name << std::endl
      << "Variables:         " << protein.potentials.size() << std::endl
      << "Edges:             " << protein.edges.size() << std::endl;
  return out;
}



/** Limit numbers to be finite */
double bound_infinity(double logvalue) {
  return logvalue;
  // // Squash the minimumlog
  // const double minLog = std::log(std::numeric_limits<double>::min())
  //   -std::log(1.0E10);
  // const double maxLog = std::log(std::numeric_limits<double>::max())
  //   -std::log(1.0E10);
  // if (logvalue < minLog) return minLog;
  // else if (logvalue > maxLog) return maxLog;
  // else return logvalue;
}





/** Parse a protein file filling in the graph g */
void parse_protein(const std::string& network_folder,
                   protein_data& protein) {
  protein.name = network_folder;
  std::string fname = network_folder + '/' + "network.bin";
  // Read in the number of variables and edges
  graphlab::binary_input_stream bis(fname.c_str());
  assert(!bis.fail());
  size_t  numvars = bis.read<int32_t>();
  size_t  numedges = bis.read<int32_t>();
  assert(!bis.fail());

  // Resize the protein data struct
  protein.potentials.resize(numvars);
  protein.edges.resize(numedges / 2);
  protein.edge_factors.resize(numedges / 2);

  
  // ignore the next 2 * num_edges * int32_t bytes this will contain
  // the list of edges which we do not need we can construct this from
  // the factors later on  
  bis.seekg(2 * numedges * sizeof(int32_t), std::ios_base::cur);  
  assert(!bis.fail());
  
  // this will now be the node potentials
  for (size_t i = 0;i < numvars; ++i) {
    size_t varid = bis.read<int32_t>();
    assert(varid == i);
    assert(!bis.fail());
    size_t cardinality = bis.read<int32_t>();
    assert(!bis.fail());
    unary_factor& factor = protein.potentials[i];
    factor.var() = i;
    factor.resize(cardinality);   
    for (size_t asg = 0; asg < cardinality; ++asg) {
      factor.logP(asg) = bound_infinity(bis.read<double>());
      assert(!bis.fail());
    }
    // Normalize the factor
    factor.normalize();
  }
  
  // this will now be the edge potentials
  size_t table_count = bis.read<int32_t>();
  assert(!bis.fail());
  assert(table_count == protein.edges.size());
  assert(table_count == protein.edge_factors.size());
  // read the binary factors
  for (size_t i = 0; i < table_count; ++i) {

    // read the src vertex
    size_t srcid = bis.read<int32_t>();
    assert(srcid < numvars);
    size_t srccard = bis.read<int32_t>();
    assert(!bis.fail());
    
    // read the destination vertex
    size_t destid = bis.read<int32_t>();
    assert(destid < numvars);
    size_t destcard = bis.read<int32_t>();
    assert(!bis.fail());

    // We only store edges in one direction
    assert(srcid < destid);
    protein.edges[i] = std::make_pair(srcid, destid);
    
    // Check that the arity match
    assert(srccard ==
           protein.potentials[srcid].arity() );
    assert(destcard ==
           protein.potentials[destid].arity() );

    // setup the vertex data blob
    protein.edge_factors[i] =
      binary_factor(srcid, srccard, destid, destcard); 
    binary_factor& factor = protein.edge_factors[i];
    for (size_t j = 0; j < srccard * destcard; ++j) {
      size_t srcasg  = bis.read<int32_t>();
      size_t destasg = bis.read<int32_t>();
      double value   = bis.read<double>();
      assert(!bis.fail());
      assert(srcasg < srccard);
      assert(destasg < destcard);
      factor.logP(srcid, srcasg, destid, destasg) =
        bound_infinity(value);
    }
    factor.normalize();
  }
}












void construct_protein_sidechain_graph(const protein_data& protein,
                                       gl::graph& graph) {
  // Initialize the vertex data to somethine sensible
  vertex_data vdata;

  // Create all the vertices 
  for(size_t i = 0; i < protein.potentials.size(); ++i) {
    // Copy the potential
    assert(i == protein.potentials[i].var());
    vdata.potential = protein.potentials[i];
    
    // vdata.asg = vdata.potential.sample();
    vdata.asg = std::rand() % vdata.potential.arity();
    
    // Initialize the belief
    vdata.belief.var() = vdata.potential.var();
    vdata.belief.resize(vdata.potential.arity());
    vdata.belief.uniform(-std::numeric_limits<double>::max());


    
    // Add the vertex
    vertex_id_t vid = graph.add_vertex(vdata);
    assert(i == vid);
  }
 
  // Add all the edges
  edge_data edata;
  for(size_t i = 0; i < protein.edges.size(); ++i) {
    size_t source = protein.edges[i].first;
    size_t target = protein.edges[i].second;
    const binary_factor& factor = protein.edge_factors[i];
    assert(source == factor.var1());
    assert(target == factor.var2());
    size_t source_arity = factor.arity1();
    size_t target_arity = factor.arity2();
    assert(source_arity == graph.vertex_data(source).potential.arity());
    assert(target_arity  == graph.vertex_data(target).potential.arity());
    assert(source == graph.vertex_data(source).potential.var());
    assert(target == graph.vertex_data(target).potential.var());
    // Initialize the edge data
    edata.weight = 1;
    edata.factor_id = i;
    // Add the forward edge
    edata.message.var() = target;
    edata.message.resize(target_arity);
    edata.message.uniform();
    graph.add_edge(source, target, edata);
    // Add the reverse edge
    edata.message.var() = source;
    edata.message.resize(source_arity);
    edata.message.uniform();
    graph.add_edge(target, source, edata);
  }
  // Finalize the graph
  graph.finalize();

  // Construct and then verify a coloring
  size_t colors = graph.compute_coloring();
  assert(graph.valid_coloring());
  std::cout << "Colors: " << colors << std::endl;
  std::vector<size_t> color_counts(colors, 0);
  for(vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) 
    color_counts[graph.color(vid)]++;
  for(size_t i = 0; i < color_counts.size(); ++i) 
    std::cout << '\t' << i << ":\t" << color_counts[i] << std::endl;
} // End of construct graph








#endif
