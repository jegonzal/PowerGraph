const RESET_PROB = 0.15;
const OUT_EDGES = 2;
const NO_EDGES = 0;

var constructor = function(){

  var last_change = 0.0;

  var gather = function(vertex, edge){
    return ((1.0 - RESET_PROB) / edge.source.num_out_edges) * edge.source.data;
  };

  var apply = function(vertex, total){
    var newval = total + RESET_PROB;
    last_change = Math.abs(newval - vertex.data);
    vertex.data = newval;
  };

  var scatter = function(context, vertex, edge){
    context.signal(edge.target);
  };

  var scatter_edges = function(vertex){
    if (last_change > 0.01) return OUT_EDGES;
    return NO_EDGES;
  };

  return {
    gather: gather,
    apply: apply,
    scatter: scatter,
    scatter_edges: scatter_edges
  };

};


var s = new pilot();
s.loadGraph("../../../demoapps/pagerank/matlab_tools/random_graph.tsv", "tsv");
s.transformVertices(function(v){ v.data = 1 });
s.fly(constructor);
s.destroy();
