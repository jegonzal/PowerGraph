#include <graphlab.hpp>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <graphlab/engine/warp_engine.hpp> 
#include <graphlab/engine/warp_parfor_all_vertices.hpp> 
#include <unsupported/Eigen/MatrixFunctions>

#include <Eigen/Dense>

int verbose;
std::vector<double> Hleft;
//#include <functional>

/*

  Total Subgraph Centrality.

  For a graph G with adjacency matrix A, TSC(G) = exp(A)*b, where 1 is
  the ones vector.


  We're going to implement this with an Arnoldi solver, following Saad
  (1992).

  The algorithm works like this:
  Choose a maximum iteration m.  Then make matrices V and H:
  b = ones(A.nodecount)
  V[0] = b/||b||
  for j in 0..m:
      w = A*V[j]
      for i in 0..j:
         H[i,j] = (w,V[i])
	 w = w - H[i,j] * V[i]
      H[j+1,j] = ||w||
      V[j+1] =w/||w||

  Then TSC = exp(A)*b ~= (V * exp(H) / ||b||)[:,0].  Stop when
  successive approximations converge, or we run out of steps.  

  We still have that matrix exponential, but it's small and dense, and
  there's an implementation in Eigen. 
  

  
*/

class node {
public: 
std::vector<double> V;
double w;
  double TSC;
  double prev;
  node(): w(0.0),TSC(0.0),prev(0.0){};
  void load(graphlab::iarchive& infile) {
    infile>>V>>w>>TSC>>prev;
  }
  void save(graphlab::oarchive& outfile) const {
    outfile<<V<<w<<TSC<<prev;
  }
  
  
};

class edge {
public:
  double weight;
  edge(): weight(1.0) {} ;
  edge(double weight) : weight(weight) {};
  

  void save(graphlab::oarchive& outfile) const {
    outfile<<weight; 
  }
  
  void load(graphlab::iarchive& infile) {
    infile>>weight; 
  }
  
  
};


typedef node vertex_data_type;
typedef edge edge_data_type;
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type; 
typedef graphlab::warp::warp_engine<graph_type> engine_type; 


// This is just a little class to be used to find the maximum change in TSC values. 
class max_finder{
public:
  double data;
  max_finder& operator+=(const max_finder& other){
    if (this->data <other.data){
      this->data = other.data;
    }
    return *this;
  }
  
  max_finder(): data(std::numeric_limits<double>::max()) {};
  max_finder(double x): data(x) {};
  void load(graphlab::iarchive& infile) {
    infile>>data;
  }
  void save(graphlab::oarchive& outfile) const {
    outfile<<data;
  }


};
  

 

// These three functions compute the A*w step of the
// Arnoldi iteration. 
double arnoldi_map(graph_type::edge_type e, 
		   graph_type::vertex_type v){
  return v.data().V.back();

}
void arnoldi_combine(double& v1, const double& v2) {
  v1 += v2;
}
void AVj_to_w(graph_type::vertex_type& v) {
  v.data().w = graphlab::warp::map_reduce_neighborhood(v, 
						       graphlab::IN_EDGES,
						       arnoldi_map,
						       arnoldi_combine);
  return;
}

// Pushes the current w onto the V matrix
void w_to_v(graph_type::vertex_type& v){
  v.data().V.push_back(v.data().w);
}
// Scales w by a constant factor.  This is meant to be called through transform_vertices 
// with a boost::bind to take care of the argument list
void scale_w(graph_type::vertex_type& v, double scale_factor){
  v.data().w /= scale_factor;
}

// Dumps the TSC to stdout
void print_TSC(graph_type::vertex_type& v){
  std::cout<<v.id()<<" " <<v.data().TSC<<std::endl;
}

// For debugging.
void print_w(graph_type::vertex_type& v){
  std::cout<<v.id()<<" " <<v.data().w<<" "<<v.data().V.size()<<" "<<v.data().V.back()<<std::endl;
}


// Make initial vector if we want a column from exp(A)
void initialize_column(graph_type::vertex_type& v, int i,int m) {
  if (v.id()==i){
    v.data().w = 1.0;
  } else{
    v.data().w = 0.0;
  }
  v.data().V.reserve(m);

}
// Make initial vector if we want the TSC
void initialize_TSC(graph_type::vertex_type& v, int m) {
  v.data().w = 1.0/sqrt((double)m);
  v.data().V.reserve(m);
}

// ||w||**2
double sum_w(const graph_type::vertex_type&  v){
  return v.data().w*v.data().w; 
}

// w*V[i].  Call via boost::bind to select the i.
double w_dot_V(const graph_type::vertex_type&  v, int i){
  return v.data().w * v.data().V[i];
}

// w - (h,w)*V[i].  Call via boost::bind to set i, hdot. 
void w_minus_hdot(graph_type::vertex_type& v, int i, double hdot) {
  v.data().w -= hdot * v.data().V[i];
}


// TSC = V*Hleft
void accumulate_hleft(graph_type::vertex_type&v){
  v.data().prev = v.data().TSC;
  v.data().TSC = 0;
  for(int j=0;j<Hleft.size();j++){
    v.data().TSC += v.data().V[j] * Hleft[j];
  }
  if (verbose){
    if (v.id()<10){
      std::cout<<"TSC "<<v.id()<<" "<<Hleft.size()<<" "<<v.data().TSC<<" "<<v.data().prev<< " "<<(v.data().TSC-v.data().prev)/(1e-15+v.data().TSC) <<std::endl;
    }
  }
  return;
}

// sum ((TSC-prevTSC)/TSC)
double total_error(const graph_type::vertex_type& v){
  return fabs((v.data().TSC-v.data().prev)/(1e-15+v.data().TSC));
}

// max ((TSC-prevTSC)/TSC)
max_finder max_error(const graph_type::vertex_type& v){
  return max_finder(fabs((v.data().TSC-v.data().prev)/(1e-15+v.data().TSC)));
}



int main(int argc, char** argv) {
  graphlab::command_line_options clopts("Total Subgraph Centrality");
  graphlab::distributed_control dc;
  std::string infile;
  std::string format = "tsv";
  int verbose=0;
  int m=100;
  int column = -1;
  clopts.attach_option("graph",infile,"Input graph.");
  clopts.attach_option("format",format,"Input format. Default tsv.");
  clopts.attach_option("m",m,"Maximum number of orthogonal vectors to approximate with. Default 100.");
  clopts.attach_option("column",column,"Column of exponential to calculate (instead of row-sum)");
  clopts.attach_option("verbose",verbose,"Verbosity level.");
  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (infile==""){
    dc.cout() <<"ERROR: Must specify --graph! \n";
    return EXIT_FAILURE;
  }
  graph_type graph(dc,clopts);
  graph.load_format(infile, format);
  graph.finalize();
  if (m>graph.num_vertices()) {
    m = graph.num_vertices();
  }
  engine_type engine(dc,graph,clopts);

  engine.signal_all();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m+1,m+1);
  Hleft.resize(m);
  for(int i=0;i<m;i++) { Hleft[i] = 0.0;}
  double beta;
  if (column>=0) {
    if (column>graph.num_vertices()){
      column = 0;
    }
    graph.transform_vertices(boost::bind(initialize_column,_1,column,m));
    beta = 1.0;
  } else {
    graph.transform_vertices(boost::bind(initialize_TSC,_1,m));
    beta = sqrt(m);
  }


  // The first column of V is just w
  graph.transform_vertices(w_to_v);


  for(int j=0;j<m;j++) {
    graphlab::warp::parfor_all_vertices(graph,AVj_to_w);
    // We should be able to move this loop inside the map-reduce and transform steps
    for(int i=0;i<=j;i++) {
      H(i,j)= graph.map_reduce_vertices<double>(boost::bind(w_dot_V,_1,i));
      graph.transform_vertices(boost::bind(w_minus_hdot,_1,i,H(i,j)));
    }
    H(j+1,j)= sqrt(graph.map_reduce_vertices<double>(sum_w));
    // If we're in this case, it means we have a spanning set and we shouldn't
    // prep for the next iteration.  Otherwise, make the data we'll need next time.
    if ( ( !std::isnan(H(j+1,j))) && (H(j+1,j) > 0) ) {
      graph.transform_vertices(boost::bind(scale_w, _1, H(j+1,j)));
      graph.transform_vertices(w_to_v); 
    }

    if (j>0) {
      // Have we converged? 
      Eigen::MatrixXd EH(H);
      EH = EH.exp();
      Hleft.resize(j+1);
      
      for(int i=0;i<=j;i++) { Hleft[i] = EH(i,0) * beta;}
      graph.transform_vertices(accumulate_hleft);
      if (j>1) {
 	max_finder largest_error= graph.map_reduce_vertices<max_finder>(max_error);
	double all_error = graph.map_reduce_vertices<double>(total_error);
	if (verbose){
	  std::cerr<<"ARNOLDI STEP FINISHED "<<j<<std::endl;
	  std::cout<<"MAX ERROR: "<<largest_error.data<<" TOTAL ERROR: "<<all_error<<std::endl;
	}
	if (largest_error.data < 1e-15) { break; };
	if (all_error < 1e-15){  break; }
      }
    }
    // If we know we are in the last iteration, break
    if (std::isnan(H(j+1,j))){ break;}
    if (fabs(H(j+1,j))<1e-15){ break;}

 
  }
  graph.transform_vertices(print_TSC);


}
