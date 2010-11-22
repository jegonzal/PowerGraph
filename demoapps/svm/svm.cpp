#include <cstdio>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>

#include <graphlab.hpp>
#include <graphlab/engine/engine_options.hpp>
#include <graphlab/macros_def.hpp>

using namespace std;
using namespace graphlab;


// Allow easy testing with float/double.
typedef double probability_t;

const size_t PRIMAL_LOSS = 0;   // double value. The estimated primal loss
const size_t NUMSVS = 1;        // integer value. The number of support vectors
const size_t TOTALGRADIENT = 2; // double. L1 gradient
const size_t STEPSIZE = 3;      // double. current stepsize
const size_t ITERATION_COUNT = 4; // double. Number of times 'b' has been updated
const size_t NUMDAT = 5;        // size_t Number of data points
const size_t NUMDIM = 6;        // size_t. Number of dimensions
const size_t BVERTEX = 7;       // size_t. vertexid of the B vertex

const double BOUND = 1E-5;
double C = 1.0;



struct dat_vertex_data {
  float primal_loss;
  char last_alpha;
  char y;
  bool sv;
};

struct feat_vertex_data {
  float primal_loss;
  float w;
  float gradient;
};

struct vertex_data {
  union {
    dat_vertex_data dat;
    feat_vertex_data feat;
  } ;
  bool isdat;
  
  void save(oarchive &arc) const{
    serialize(arc, this, sizeof(vertex_data));
  }
  void load(iarchive &arc) {
    deserialize(arc, this, sizeof(vertex_data));
  }
};


struct edge_data {
  float x;
};


double get_loss(const vertex_data &v) {
  if (v.isdat) {
    return v.dat.primal_loss;
  }
  else {
    return v.feat.primal_loss;
  }
}

size_t get_sv(const vertex_data &v) {
  if (v.isdat) return size_t(v.dat.sv);
  else return 0;
}

double get_gradient(const vertex_data &v) {
  if (v.isdat) return 0;
  else return fabs(v.feat.gradient); 
}

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;

void dat_update_function(gl_types::iscope &scope,
                         gl_types::icallback &scheduler,
                         gl_types::ishared_data* shared_data);
void w_update_function(gl_types::iscope &scope,
                       gl_types::icallback &scheduler,
                       gl_types::ishared_data* shared_data);

                       
bool primal_terminator(const gl_types::ishared_data* shared_data) {
  assert(shared_data != NULL);
  double stepsize = shared_data->get(STEPSIZE).as<double>();
  double l1gradient = shared_data->get(TOTALGRADIENT).as<double>();
  size_t numdim = shared_data->get(NUMDAT).as<size_t>();
  std::cout << stepsize*(l1gradient/numdim) << "\n";
  return stepsize*(l1gradient/numdim) <= BOUND;
}


/// ====== UPDATE FUNCTION ======= ///
void dat_update_function(gl_types::iscope &scope,
                         gl_types::icallback &scheduler,
                         gl_types::ishared_data* shared_data) {
  
  size_t numdat = shared_data->get(NUMDAT).as<size_t>();
  dat_vertex_data& thisvertex = scope.vertex_data().dat;
  // compute alpha
  
  double s = 0;
  foreach(edge_id_t eid, scope.in_edge_ids()) {        
    const edge_data& edata = scope.edge_data(eid);  
    const feat_vertex_data& feat = scope.neighbor_vertex_data(scope.source(eid)).feat;
    s = s + edata.x * feat.w;
  }
  s = s * thisvertex.y;

  thisvertex.primal_loss = (C/numdat * std::max(0.0, 1-s));

  char alpha = -1;
  if (s>=1) {
    alpha = 0;
    thisvertex.sv = 0;
  }
  else thisvertex.sv = 1;

  if (thisvertex.last_alpha != alpha) {
    thisvertex.last_alpha = alpha;
    // schedule all connected w's if I change
    foreach(edge_id_t eid, scope.in_edge_ids()) {
      const edge_data& edata = scope.edge_data(eid);
      scheduler.add_task(gl_types::update_task(scope.source(eid),
                                               w_update_function), edata.x);
    }
  } else if (alpha != 0) {
    // schedule all connected w's if alpha is non-zero
    foreach(edge_id_t eid, scope.in_edge_ids()) {
      const edge_data& edata = scope.edge_data(eid);  
      scheduler.add_task(gl_types::update_task(scope.source(eid),
                                               w_update_function), edata.x);
    }
  } 
}

void w_update_function(gl_types::iscope &scope,
                       gl_types::icallback &scheduler,
                       gl_types::ishared_data* shared_data) {
  // compute the local gradient
  // am I b?
  feat_vertex_data& thisvertex = scope.vertex_data().feat;
  double gradient = 0.0;
  double stepsize = shared_data->get(STEPSIZE).as<double>();
  size_t numdat = shared_data->get(NUMDAT).as<size_t>();
  size_t bvertex = shared_data->get(BVERTEX).as<size_t>();
  
  if (scope.vertex() != bvertex) {

    gradient = thisvertex.w;
    foreach(edge_id_t eid, scope.out_edge_ids()) {        
      const edge_data &edata = scope.edge_data(eid);
      const dat_vertex_data& dat = 
            scope.neighbor_vertex_data(scope.target(eid)).dat;
      if (dat.last_alpha != 0) {
        gradient = gradient + C * dat.y * dat.last_alpha*edata.x;
      }
    }
    // take a gradient step
    gradient = -gradient;
    gradient = gradient / numdat;
  } else {
    foreach(edge_id_t eid, scope.out_edge_ids()) {        
      const dat_vertex_data& dat =
            scope.neighbor_vertex_data(scope.target(eid)).dat;
      if (dat.last_alpha != 0) {
        gradient = gradient - C * dat.y * dat.last_alpha;
      }
    }
    // take a gradient step
    gradient = -gradient;
    gradient = gradient / numdat;
  }
  
  if (scope.vertex() == bvertex) {
    shared_data->atomic_apply(ITERATION_COUNT, 
                              gl_types::apply_ops::increment<size_t>, 
                              size_t(1));
    size_t iterationcount = shared_data->get(ITERATION_COUNT).as<size_t>();
    stepsize = 100.0 / (C*iterationcount);
    shared_data->atomic_set(STEPSIZE, stepsize);
/*    shared_data->sync(NUMSVS);
    shared_data->sync(PRIMAL_LOSS);
    shared_data->sync(TOTALGRADIENT);*/
  }
  thisvertex.primal_loss = thisvertex.w * thisvertex.w / 2;
  thisvertex.gradient = gradient;
  //std::cout << scope.vertex() << ": " << thisvertex->w << "\n";
  if (std::fabs(stepsize * gradient) >= BOUND) {
    thisvertex.w += stepsize * gradient;
    // schedule all vertices which might change
    foreach(edge_id_t eid, scope.out_edge_ids()) {
      const edge_data& edata = scope.edge_data(eid);
      const dat_vertex_data& dat =
        scope.neighbor_vertex_data(scope.target(eid)).dat;
      if (dat.last_alpha == 0) {
        // Is the gradient in the right direction to potentially push alpha to -1?
        if (gradient * edata.x * dat.y < 0) {
          // oops shift in the wrong direction
          // schedule it
          scheduler.add_task(gl_types::update_task(scope.target(eid),
                                                   dat_update_function),
                             fabs(gradient * edata.x));
        }
      } else if (dat.last_alpha == -1) {
        // Is the gradient in the right direction to potentially push alpha to -1?
        if (gradient * edata.x * dat.y > 0) {
          // oops shift in the wrong direction
          // schedule it
          scheduler.add_task(gl_types::update_task(scope.target(eid),
                                                   dat_update_function),
                             fabs(gradient * edata.x));
        }
      }
    }
    // reschedule myself
    scheduler.add_task(gl_types::update_task(scope.vertex(), w_update_function),
                       std::fabs(gradient));
  }
}

void create_graph(graph_type& g,
                  gl_types::thread_shared_data &sdm,
                  std::string infile) {
  std::ifstream fin(infile.c_str());
  // first pass. Count the number of data points
  size_t numdat = 0;
  size_t numdim = 0;
  size_t numedges = 0;
  while (fin.good()) {
    size_t dat,feat;
    double val;
    fin >> dat>>feat>>val;
    numdat = std::max(numdat,dat);
    numdim= std::max(numdim,feat);
    ++numedges;
  }
  numdat++;
  // entry 0 is class
  fin.close();
  // add the feature vertices
  for (size_t i = 0;i < numdat; ++i) {
    vertex_data d;
    d.dat.last_alpha = 0 ;
    d.dat.y = 0;
    d.dat.primal_loss = 1000;
    d.dat.sv = true;
    d.isdat = true;
    g.add_vertex(d);
  }
  std::cout << numdim << " dimensions\n";
  std::cout << numdat << " data points\n";
  for (size_t i = 0;i < numdim; ++i) {
    vertex_data d;
    d.feat.w = 0 ;
    d.feat.primal_loss = 1000;
    d.feat.gradient = 1000;
    d.isdat = false;    
    g.add_vertex(d);
  }
  // add the 'b' vertex
  vertex_data b;
  b.feat.w = 0 ;
  b.isdat = false;
  size_t bvertex = g.add_vertex(b);

  // add an edge from bvertex to all the data
  for (size_t i = 0;i < numdat; ++i) {
    edge_data e;
    e.x=-1.0; // allows the datavertex to 'ignore' the fact that
    // this edge is special in the alpha computation
    g.add_edge(bvertex, i, e);
  }
  
  size_t numread = 0;
  fin.open(infile.c_str()); 
  while (fin.good()) {
    size_t dat,feat;
    double val;
    fin >> dat>>feat>>val;
    if (fin.fail()) break;
    numread++;
    if (numread % 100000==0) {
      std::cout <<numread << "/" << numedges << "\n";
    }
    // this is class
    if (feat == 0) {
      if (!(val ==-1 || val == 1)) {
        std::cout << "Classes should be -1 or 1 \n";
        std::cout << dat << " " << feat << " " << val << "\n";
        val = -1;
      }
      g.vertex_data(dat).dat.y = val;      
    }
    else {
      // add an edge from feature to the data
      feat--; // remember to drop the class feature
      edge_data e;
      e.x=val;
      g.add_edge(feat+numdat, dat, e);
    }    
  }
  std::cout << "Graph Read\n";
  fin.close();
  g.finalize();

  
  // initialize the shared data
  sdm.set_sync(PRIMAL_LOSS, 
               gl_types::sync_ops::sum<double, get_loss>,
               gl_types::apply_ops::identity<double>,
               double(0),
               100);

  sdm.set_sync(NUMSVS, 
               gl_types::sync_ops::sum<size_t, get_sv>,
               gl_types::apply_ops::identity<size_t>,
               size_t(0),
               100);
  
  sdm.set_sync(TOTALGRADIENT, 
               gl_types::sync_ops::sum<double, get_gradient>,
               gl_types::apply_ops::identity<double>,
               double(0),
               100);
               
  sdm.create_atomic(STEPSIZE,  double(1.0));
  sdm.create_atomic(ITERATION_COUNT, size_t(0));
  sdm.create_atomic(NUMDAT, numdat);
  sdm.create_atomic(NUMDIM, numdim);
  sdm.create_atomic(BVERTEX, bvertex);

  std::cout << g.num_vertices() << " vertices\n";
  std::cout << g.num_edges() << " edges\n";
}



bool classify(graph_type &g, vertex_id_t vertexid){
  dat_vertex_data& thisvertex = g.vertex_data(vertexid).dat;
  // compute alpha
  double s = 0;
  foreach(edge_id_t eid, g.in_edge_ids(vertexid)) {        
    edge_data edata = g.edge_data(eid);
    feat_vertex_data& feat = g.vertex_data(g.source(eid)).feat;
    s = s + edata.x * feat.w;
  }
  if (s>0 && thisvertex.y>0) return true;
  if (s<0 && thisvertex.y<0) return true;
  return false;
}

void start(gl_types::iengine *engine, graph_type &g, 
           gl_types::thread_shared_data &sdm) {
  
  size_t numdat = sdm.get(NUMDAT).as<size_t>();
  // schedule all the data vertices
  engine->add_terminator(primal_terminator);
  for (size_t i = 0;i < numdat; ++i) {
    engine->get_scheduler().
      add_task(gl_types::update_task(i, dat_update_function), 1000.0);
  }
  timer ti;
  ti.start();
  engine->start();

  
  size_t correctcount = 0;
  for (size_t i = 0;i < numdat; ++i) {
    correctcount += classify(g, i);
  }
  std::cout << "w = ";
  for (size_t i = numdat ;i < g.num_vertices(); ++i) {
    std::cout << g.vertex_data(i).feat.w << " ";
  }
  std::cout << "\n";
  sdm.sync(PRIMAL_LOSS);
  std::cout << correctcount << "/" << numdat << "\n";;
  std::cout << "Primal Loss: " << sdm.get(PRIMAL_LOSS).as<double>() << "\n";
  std::cout << "NumSVS : " << sdm.get(NUMSVS).as<size_t>() << "\n";
  std::cout << "Gradient: " << sdm.get(TOTALGRADIENT).as<double>() << "\n";
  std::cout << ti.current_time() << " seconds\n";
}

  
  

int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);


  std::string infile;
  std::string outfile;

  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("SVM input file");
  // Set the program options
  desc.add_options()
  ("infile",  boost_po::value<std::string>(&(infile))->default_value(""),
   "Input filename ")
   ("outfile",  boost_po::value<std::string>(&(outfile))->default_value(""));

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);

  if (infile == "") {
    std::cout << "Input Graph needed\n";
    return 0;
  }
  graph_type graph;
  gl_types::thread_shared_data sdm;
  create_graph(graph, sdm, infile);
  
  command_line_options clopts;
  clopts.parse(argc, argv);
  gl_types::iengine *engine = clopts.create_engine(graph);
  engine->set_shared_data_manager(&sdm);
  sdm.sync(PRIMAL_LOSS);
  sdm.sync(NUMSVS);
  sdm.sync(TOTALGRADIENT);

  start(engine, graph, sdm);
  
  
}

