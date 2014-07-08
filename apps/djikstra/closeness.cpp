
class ClosenessAlgorithm :
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE {
    bool changed;

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
	Gather g;
	for(iter*){
		long key=0;
		if((edge.source().launched == true)&&(edge.source().done == false)){
			double c = edge_type.data() + edge.source().data().cost;
			g.cost = c;
			g.id = edge.source().data().id;
		}else{
			g.id=0;
		}
		return *g;	
	}
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
	for(iter *){
		long key = 0;
		if(vertex.data()[key].launched = false){
			vertex.data()[key].launched = true;
      			if(vertex.data()[key].cost > total.data().cost){
				vertex.data()[key].cost = total.data().cost;
				vertex.data()[key].id = total.data().id;
      			}else{
				vertex.data()[key].done = true;
			}
      		}else{
			vertex.data()[key].done = true;
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
                }
        }
  };


