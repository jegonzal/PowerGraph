const RESET_PROB = 0.15;
const OUT_EDGES = 2;
const NO_EDGES = 0;

var PageRank = (function(){

  var cls = function(){

    var last_change = 0.0;

    this.gather = function(vertex, edge){
      return ((1.0 - RESET_PROB) / edge.source.num_out_edges) * edge.source.data;
    };

    this.apply = function(vertex, total){
      var newval = total + RESET_PROB;
      last_change = Math.abs(newval - vertex.data);
      vertex.data = newval;
    };

    this.scatter = function(context, vertex, edge){
      context.signal(edge.target);
    };

    this.scatter_edges = function(vertex){
      if (last_change > 0.01) return OUT_EDGES;
      return NO_EDGES;
    };

  };

  return cls;

})();


var s = new pilot();
s.loadSyntheticPowerlaw(1000);
s.transformVertices(function(v){ v.data = 1.0; });
s.fly(PageRank);
s.destroy();
