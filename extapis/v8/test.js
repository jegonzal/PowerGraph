const RESET_PROB = 0.15;

var PageRank = (function(){

  var cls = function(){

    var last_change = 0.0;

    this.gather = function(vertex, edge){
      return ((1.0 - RESET_PROB) / edge.source().numOutEdges()) * edge.source().data;
    };

    this.apply = function(vertex, total){
      var newval = total + RESET_PROB;
      last_change = Math.abs(newval - vertex.data);
      vertex.data = newval;
    };

    this.scatter = function(context, vertex, edge){
      context.signal(edge.target());
    };

    this.scatter_edges = function(vertex){
      if (last_change > 0.01) return OUT_EDGES;
      return NO_EDGES;
    };

  };

  return cls;

})();


var s = new pilot();
s.loadGraph("../../../demoapps/pagerank/matlab_tools/random_graph.tsv", "tsv");
s.transformVertices(function(v){ v.data = 1.0; });
s.fly(PageRank);
s.saveGraph("out.tsv",
  function(v){
    return v.id() + "\t" + v.data + "\n";
  },
  null, false, true, false);
s.destroy();
