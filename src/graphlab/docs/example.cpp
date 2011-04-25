/**
\page pagerank_example PageRank Example

\section sec_pagerank_example_preliminaries Preliminaries

First we need to create C++ file for our application and configure it for compilation. This is easy.

\b 1. Create directory pagerank under directory apps: cd apps mkdir pagerank 

\b 2. Create file CMakeLists.txt in the new directory and add following lines to it:
\verbatim
project(GraphLab) 
add_executable(pagerank pagerank.cpp)   
\endverbatim

\b 3. Then create file pagerank.cpp with following content: 
\code
/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */

#include <string>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");
   
  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm.");

  // Create a graphlab core
  gl_types::core core;
  
  // Parse arguments
  if(!clopts.parse(argc, argv)) {
     std::cout << "Error in parsing input." << std::endl;
     return EXIT_FAILURE;
  }
  
  // Set the engine options
  core.set_engine_options(clopts);
  
  // Create a synthetic graph
  //create_graph(core.graph());

  // Schedule all vertices to run pagerank update on the
  // first round.
  // core.add_task_to_all(pagerank_update, 100.0);
  
  // Run the engine
  double runtime = core.start();
  
  return EXIT_SUCCESS;
} 
\endcode



Above is just the bare skeleton for the application, and contains 
the code for initializing graphlab. We will soon add more code 
to include the actual functionality.

\b 4. Then finally, we need to finish the setup and check we 
can compile our new app. Go to the root of graphlab and run:

\verbatim
./configure
\endverbatim

Then, to compile the new app go to directory debug/apps/pagerank 
and run \c make . You should compilation errors, but do not 
worry about them now! 


\section sec_pagerank_example_algo Short introduction to PageRank algorithm

The objective of PageRank algorithm is to assign a score of the "importance" of
each webpage (or scientific article).
The score (Pagerank) of each page depends on the pageranks of pages <em>linking to it</em>.
In essence, the PageRank of a page is a weighted sum of pageranks of articles linking to it.  

That is, each page gives 1/N of its pagerank to each page its links to,
where N is the number of links from the page. This can be interpreted by thinking
that there is probablity of 1/N to move to any of the linked pages. 
In addition, for a more realistic ranking, a certain weight, depending on a <em>damping factor</em> is given to all
other pages in the Web. This models the situation when surfer moves to a page
not linked from a page, for example by selecting a bookmark.

In this Tutorial, we implement PageRank inference algorithm 
using a <em>power iteration method</em>. At each iteration, 
we set the pagerank of each page to be 
the weighted average of pageranks of pages linking to the page, 
plus the fraction of possibility to moving to a random page. 
Algorithm is terminated when
the pageranks change less than a predefined threshold. The exact algorithm will
become clear when we implement the <b>update function</b>.

For more information, please consult  the Wikipedia article about 
PageRank ( http://en.wikipedia.org/wiki/PageRank )



\section sec_pagerank_example_datagraph Creating the data graph

The data model of this algorithm is readily in the form of a graph.
Each page is represented by a <b>vertex</b>, of which data contains simply
the current estimate of the rank. In addition, since GraphLab graphs do not
support self-edges, we store the weight of the possible self-edge as a field
in the vertex. 

Each page has an outgoing <b>edge</b> to each vertex representing a page it has a link to. 
<b>Edge data</b> consists of a weight-attribute, which is 1/N for each outgoing edge 
(N was the number of outgoing edges of a vertex).

In addition, each edge contains a field for storing the <em>previous pagerank</em> used
from the <em>neighboring</em> vertex. This is used for computing the change of values 
between iterations, and is used to determined when to terminate. As an example, 
let there be a link between vertex 1 and 3. At time T, vertex 3 is updated and it 
reads the pagerank of vertex 1, which is <tt>0.5</tt>. It writes
this value to the edge from 1 to 3. Then, at timestep T+k, vertex 1 is updated
and it computes new pagerank <tt>0.6</tt>. Since <tt>|0.6-0.5| = 0.1</tt> is 
larger than our predefined threshold (0.00001 for example), the update function
asks the scheduler to update vertex 3 again. 

Ok, we are now ready to define our graph data types and create a simple graph.
Later on, we give also code and data for a larger graph to be loaded from a 
file. 

<b>1.</b> First we define structs for our vertex and edge datatypes. 
Add following code before the <tt> main()</tt> Function:

\code
//
// Stores the pagerank and the self weight
//
struct vertex_data {
  float value;
  float self_weight; 
  vertex_data(float value = 1) : value(value), self_weight(0) { }
}; 


//
// Edge data represents the weight as well as the weight times the
// last value of the source vertex when the target value was computed.
//
struct edge_data {
  float weight;
  float old_source_value;
  edge_data(float weight) :
    weight(weight), old_source_value(0) { } 
    
}; 
\endcode


<b>2.</b> Since GraphLab uses C++ templates for convenience and for 
avoiding casting errors, we need to tell the compiler the type of  our graph, which we
name <tt>pagerank_graph</tt>. In addition, we declare that all GraphLab classes
will be templated based on the <tt>pagerank_graph</tt>. Add following code after
the previous definitions of vertex and edge data structures.

\code
// The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> pagerank_graph;

//
// The collection of graphlab types restricted to the graph type used
// in this program.
//
typedef graphlab::types<pagerank_graph> gl_types;
\endcode



<b>3.</b> To introduce basics of graph creation, we construct a 5-vertex (page) 
graph by hand.  In the following code snippet, function <tt>create_graph()</tt> first creates
five vertices with default vertex data. Vertices will get ids from 0 to 4. After
that, a set of edges are created. In addition, the fifth page (id 4) links to
all pages, including itself, and thus needs to set the <tt>self_weight</tt> 
parameter of the vertex data. Please note that each edge has a weight that is
inverse of the number of edges from the vertex it is from.

\code
// Creates simple 5 vertex graph
void create_graph(pagerank_graph& graph) {
	// Create 5 vertices
	graph.add_vertex(vertex_data());
	graph.add_vertex(vertex_data());
	graph.add_vertex(vertex_data());
	graph.add_vertex(vertex_data());
	graph.add_vertex(vertex_data());
	
	// Page 0 links to page 3 only, so weight is 1
	graph.add_edge(0, 3, edge_data(1));
	
	// Page 1 links to 0 and 2
	graph.add_edge(1, 0, edge_data(0.5));
	graph.add_edge(1, 2, edge_data(0.5));
	
	// ... and so on
	graph.add_edge(2, 0, edge_data(1.0/3));
	graph.add_edge(2, 1, edge_data(1.0/3));
	graph.add_edge(2, 3, edge_data(1.0/3));
	graph.add_edge(3, 0, edge_data(0.25));
	graph.add_edge(3, 1, edge_data(0.25));
	graph.add_edge(3, 2, edge_data(0.25));
	graph.add_edge(3, 4, edge_data(0.25));
  graph.add_edge(4, 0, edge_data(0.2));
 	graph.add_edge(4, 1, edge_data(0.2));
	graph.add_edge(4, 2, edge_data(0.2));
	graph.add_edge(4, 3, edge_data(0.2));
	// and self edge from 4 to 4, which must be handled specially.
	graph.vertex_data(4).self_weight = 0.2;
}
\endcode


<b>4.</b> Then <i>uncomment</i> following line from the <tt>main()</tt>:


\code
  // create_graph(core.graph());
\endcode


Graphlab object <tt>core</tt> contains handle to the <b>graph</b>, computation
<b>engine</b> and allows configuration of the runtime. Please see section 
\ref intro_sec for more information.


\section sec_pagerank_example_update_function Update function


Now we are ready to write the actual computation, which was described above.
To accomplish this, we define an <b>update function</b> <tt>pagerank_update</tt>. 
Update function operates on one vertex a time. A typical GraphLab update function consists
of three steps:

\li 1. Read values from (inbound) neighbor vertices and edges; accumulate 
       into a local variable to get a new value.
\li 2. Write the new value to the vertex data.
\li 3. If new value changed more than a threshold, ask scheduler to update the
(outbound) neighbors vertices. Optionally set a priority for the update. 

Algorithm will end when no update function decides to add neighbors
again. It is important to note, that if a vertex was already scheduled by 
another neighbor, it is not added again to the task list (we call this
<i>task pruning</i>).

<b>1.</b> First we define two constant parameters for the algorithm: the damping
factor and termination threshold. Usually, it is better to use  
<b>shared variables</b> for managing such paramaters, but to keep this Tutorial simple,
we "hard-code" them.  See the supplementary sections below how to use the SDT for this purpose.

Add following lines somewhere in the beginning of the file.

\code
#define termination_bound 1e-5
#define damping_factor 0.85   // PageRank damping factor
\endcode


<b>2.</b> Now let's add the empty update function, before the <tt>main()</tt> function:

\code
//
// The PageRank update function
//
void pagerank_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler) {
  
  }
\endcode
 

All update functions share the same signature. Note that the parameters are
passed as <i>references</i>.



<b>3.</b> Now, we start by adding functionality to the update function. First we
assign a reference of the vertex the scope is operating on to variable <tt>vdata</tt>:


\code
   // Get the data associated with the vertex
  vertex_data& vdata = scope.vertex_data();
\endcode


<b>4.</b> Then we add the main computation loop of the update function. 
We loop over each <i>incoming</i> edge, and read the pagerank value of the 
associated vertex and multiply it by the weight of the edge. These values are
accumulated to local variable <tt>sum</tt>. Immediatelly after reading the value, we record
the read value to the edge. This is later used, when the neighboring vertex is
updated, to determine the change between previous used value.


\code
  
  // Start with including contribution of the self-edge
  float sum = vdata.value * vdata.self_weight;

  foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
    // Get the neighobr vertex value
    const vertex_data& neighbor_vdata =
      scope.const_neighbor_vertex_data(scope.source(eid));
    double neighbor_value = neighbor_vdata.value;
    
    // Get the edge data for the neighbor
    edge_data& edata = scope.edge_data(eid);
    // Compute the contribution of the neighbor
    double contribution = edata.weight * neighbor_value;
    
    // Add the contribution to the sum
    sum += contribution;
    
    // Remember this value as last read from the neighbor
    edata.old_source_value = neighbor_value;
  }
\endcode


<b>5.</b> Now we are ready to compute the new PageRank for the vertex 
(page), which is a weighted sum of the computed <tt>sum</tt> and the 
<i>damping factor</i>. New value is stored in field <tt>value</tt> of the vertex.

\code
  // compute the new pagerank
  sum = (1-damping_factor)/scope.num_vertices() + damping_factor*sum;
  vdata.value = sum;
\endcode


<b>6.</b> Now we need to check if the new value differed from previous value
significantly. If it did, we ask the scheduler to update the neighbor vertices.
Note that we check the change for each neighbor separately, by computing
a weighted <i>residual</i> (this is the most
efficient way, but requires a lot of memory if the graph is big, because the previous value needs
to be stored in each edge separately. Simpler, but less accurate way, is to
just compare old and new value of the vertex). 

\code
// Schedule the neighbors as needed
  foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
    edge_data& outedgedata = scope.edge_data(eid);
    
    // Compute edge-specific residual by comparing the new value of this
    // vertex to the previous value seen by the neighbor vertex.
    double residual =
      outedgedata.weight *
      std::fabs(outedgedata.old_source_value - vdata.value);
    // If the neighbor changed sufficiently add to scheduler.
    if(residual > termination_bound) {
      gl_types::update_task task(scope.target(eid), pagerank_update);
      scheduler.add_task(task, residual);
    }
  }
\endcode


<b>7.</b> Finally, we need to kick-start the computation by scheduling each
vertex for update. To achieve this, <i>uncomment</i> following line from the <tt>main()</tt> function:

\code
   // core.add_task_to_all(pagerank_update, 100.0);
\endcode


\section sec_pagerank_example_output Output


After the algorithm has converged, we would like to output the results. 
In addition, to get comparable results, we should normalize the results.
The best way would be to use <b>sync</b> facility, but for simplicity 
we do it manually.



<b>1.</b> To compute the normalizer (i.e sum of pageranks), please add following
code in the end of the <tt>main()</tt> function, before <tt>return</tt>:


\code
  double norm = 0.0;
  for(graphlab::vertex_id_t vid=0; vid&lt;core.graph().num_vertices(); vid++) {
  	 norm += core.graph().vertex_data(vid).value;
  }
\endcode


<b>2.</b> Then we simply output normalized pagerank of five first vertices:


\code
 for(graphlab::vertex_id_t vid=0; vid&lt;5 ; vid++) {
  	 std::cout << "Page " << vid << " pagerank = " <<
  	 	core.graph().vertex_data(vid).value / norm << std::endl;
  }
\endcode

\section sec_pagerank_example_running_code Running the code

<b>1.</b> To compile, do as before and go to <tt>debug/apps/pagerank</tt> and
enter <tt>make pagerank</tt>. For better code performance, you can instead
go to <tt>release/apps/pagerank</tt>.

<b>2.</b> Now you can run the algorithm by simply entering <tt>./pagerank</tt>.

You should see approximately following output:

\verbatim
Graphlab finished, runtime: 0.000952 seconds.
Page 0 pagerank = 0.235752
Page 1 pagerank = 0.165445
Page 2 pagerank = 0.183704
Page 3 pagerank = 0.301708
Page 4 pagerank = 0.11339
\endverbatim



<b>3.</b> You can also try running with different schedulers, number of threads and consistency
settings. With such a small data, you should not expect to see any difference
in performance, but it is good to learn how to declare runtime parameters. 
Here are some examples:

\verbatim
./pagerank --scheduler=multiqueue_fifo
./pagerank --ncpus=1
./pagerank --ncpus=8
./pagerank --ncpus=16
./pagerank --scope=edge
./pagerank --scope=full
\endverbatim




\section sec_pagerank_example_sup1 Supplementary: loading  graphs from file

Following code reads a dataset for PageRank from a file. The format
is explained in the header. To use it, simply replace the
call to <tt>create_graph()</tt> in the main with
<tt> load_graph("mydata.txt", core.graph());</tt>

\code
/*
 * Load a graph file specified in the format:
 *
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph(const std::string& filename,
                pagerank_graph& graph) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    fin.ignore(1); // skip comma
    fin >> target;
    fin.ignore(1); // skip comma
    fin >> weight;
    //    fin.ignore(1); // skip comma
    // Ensure that the number of vertices is correct
    if(source >= graph.num_vertices() ||
       target >= graph.num_vertices())
      graph.resize(std::max(source, target) + 1);
    if(source != target) {
      // Add the edge
      edge_data edata(weight);
      graph.add_edge(source, target, weight);
    } else {
      // add the self edge by updating the vertex weight
      graph.vertex_data(source).self_weight = weight;
    }       
  }
  graph.finalize();
  return true;
} 
\endcode


\section sec_pagerank_example_sup2 Supplementary: using shared variables for parameters

Instead of defining termination threshold and damping factor as
hard-coded constants, it is better to manage them with the
<b>Shared Variables</b>. This makes it easy to pass them from
command-line, for example. Following code achieves this:

First define the variables:
\code
// Keys for shared data 
glshared_const<float> DAMPING;
glshared_const<float> TERMINATION_BOUND;
\endcode


Then before starting graphlab, assign values for the keys as constants:
\code
DAMPING.set(0.85);
\endcode


Then in the update function, you can access the values:
\code
 float damping_factor = DAMPING.get_val();
\endcode

\section sec_pagerank_example_sup3 Supplementary: adding command line parameters

GraphLab uses Boost Program Options package for command line
option parsing. You can add your own options, such as <tt>--dampingfactor=xx</tt>
easily. Following code shows an example how to add one custom string parameter for input file
and float parameter for the termination bound:


\code
  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm on a input file.");
  clopts.attach_option("infile", &filename,
                       "PageRank input file. In src, dest, weight format.");
  clopts.attach_option("bound", &termination_bound, termination_bound,
                       "Termination bound for maximum change in vertex values.");
  clopts.add_positional("infile");
  
  if(!clopts.parse(argc, argv)) {
     std::cout << "Error in parsing input." << std::endl;
     return EXIT_FAILURE;
  }
    
  // Set the engine options
  core.set_engine_options(clopts);
\endcode


*/
