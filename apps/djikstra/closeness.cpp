
class BetweenessGather{
public:
	map<long,long> counts;
	map<long,long> edge_count;
	
 void save(graphlab::oarchive& oarc) const {
    oarc << counts << active;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> counts >> active;
  }

  BetweenessGather& operator+=(const& BetwenesisGather gather){
	for(iter*){
                long key = 0;
        	this.counts[key] += gather.counts[key];
		this.edge_count[key] += gather.edge_count[key];
	}
        for(iter*){
                long key = 0;
                if(this.counts.contains(key)){
                        this.counts[key] = gather.counts[key];
			this.edge_count[key] = gather.edge_count[key];
                }
        }
	return this;
	
};

typedef BetweenessGather gather_type;

class BetweennessAlgorithm :
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
	BetweenessGather g;
	for(iter*){
		long key=0;
		if(edge.source().launched == true){
			g.map[key] = edge.source().data().map[key].count;
			g.count = 1;
		}
	}
	return *g;	
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
	for(iter *){
		long key = 0;
		if((vertex.data()[key].launched == true)&&(vertex.data()[key].done == false)&&(vertex.data()[key].count==total.edge_count[key])){
			vertex.data()[key].launched = true;
			vertex.data()[key].count = total[key].count;
      		}
	}
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      // if vertex data changes, scatter to all edges.
     	bool done = true;
	for(iter* ){
		long key = 0;
		if(vertex_type.data()[key].launched && !vertex_type.data()[key].done){ 
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
        for(iter* ){
                long key = 0;
                if((vertex.data()[key].done == false) && (vertex.data()[key].launched == true)){
                        context.signal(edge.target());
			vertex.data().map[key].done = true;
                }
        }
  };


