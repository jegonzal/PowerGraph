/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
\page detailed_example Detailed Example

In this section we provide a brief tutorial on how the pieces of 
GraphLab come together to form a simple GraphLab program. In this 
tutorial we construct a synthetic application which uses many of 
the GraphLab concepts. However, many real GraphLab applications 
will not need all the pieces described in this tutorial. For a 
syntax highlighted (condensed version) of the demo.cpp see demo.cpp

\section detailed_example_problem The Problem

We begin by first describing the synthetic problem we will be solving. We are given a 
d by d grid of magnetic coins. Each coin only interacts with its 4 neighboring coins 
as illustrated by the edges in Figure 1(a). When a coin is flipped one of the following 
two outcomes occurs:

\li the coin lands red with probability proportionate to the number of red neighbors
\li the coin lands black with probability proportionate to the number of black neighbors

\htmlonly
  <center>
    <table width="600px">
      <tr> 
        <td align="center">
          <img src="fig_02.gif" alt="Mixed Grid" width="160"/>
        </td>
        <td align="center">
          <img src="fig_03.gif" alt="Mixed Grid" width="160"/>
        </td>
        <td align="center">
          <img src="fig_04.gif" alt="Mixed Grid" width="160"/>
        </td>
      </tr>
      <tr>
        <td align="center"><b>(a)</b></td> <td align="center"><b>(b)</b></td> <td align="center"><b>(c)</b></td>
      </tr>
      <tr>
      <td colspan="3">
      This picture illustrates the grid of red/black coins.  Each coin
      is connected to its 4 cardinal neighbors (as shown by the
      edges). <b>(a)</b> Illustrates an intermediate state of the system. 
      <b>(b,c)</b> The two stable states of the system.
      </td>
      </tr>
    </table>
  </center>
\endhtmlonly

The only two stable joint states are when all coins are black, Figure 1(b), 
or where all coins are red, Figure 1(c). We are interested in computing the 
average number of flips for each coin before a stable state is reached. 

\section detailed_example_includes Includes

We first need to include some headers: <tt>graphlab.hpp</tt> includes everything 
you will ever need from GraphLab. 

\code
// standard C++ headers 
#include <iostream> 
// includes the entire graphlab framework 
#include <graphlab.hpp> 
\endcode

In some of our demo code you may occasionally encounter the following additional included header: 

\code
#include <graphlab/macros_def.hpp> 
\endcode

If you use <tt>macros_def.hpp</tt> in a header file, it must be paired with a matching 
<tt> \#include <graphlab/macros_undef.hpp> </tt> at the end of the header file.

\section detailed_example_data_graph Getting Started with The Data Graph 

The first step to designing a GraphLab program is setting up the data graph. 
To do this we will need to define the data elements and their dependencies. 
The primary data element in this simple program is the data at each vertex 
which records the number of flips so far and the current color of the vertex. 
Here we will assume the color red is true and the color black is false. 

\code
struct vertex_data { 
  size_t numflips; 
  bool color; 
}; 
\endcode

GraphLab provides facilities to directly save/load graphs from disk. 
However, to do so,  you must implement a save and load function
so that GraphLab can understand your data-structures.  
The serialization mechanism is simple to use and it understands all basic datatypes 
as well as standard STL containers. If the STL container contains non-basic datatypes 
(such as a struct), save/load functions must be written for the datatype. If we wanted 
to be able to save the graph to a file we would implement the following functions.
See \ref Serialization for more details.
\code
struct vertex_data { 
  size_t numflips; 
  bool color; 

  void save(graphlab::oarchive& archive) const { 
    archive << numflips << color; 
  } 

  void load(graphlab::iarchive& archive) { 
    archive >> numflips >> color; 
  } 
};
\endcode

In this example, we do not need edge data. However GraphLab currently 
does not have a mechanism to completely disable the use of edge data. 
Therefore, we will just put an arbitrary small placeholder type on the edges.
\code
typedef char edge_data;
\endcode
Note that we do not need to write a save/load function here since 
GraphLab's serializer already understands basic datatypes. One could imagine
an alternative problem where edge weights are associated with each pair of 
interacting coins. 


\section detailed_example_typedefs GraphLab Typedefs

The GraphLab graph is templatized over the vertex data as well as the edge data. 
Here we define the type of the graph using a typedef for convenience.
\code
typedef graphlab::graph<vertex_data, edge_data> graph_type;
\endcode
Since graphlab is heavily templatized and can be inconvenient to use in its standard form, 
the <tt>graphlab::types</tt> structure provides convenient typedefed "shortcuts" to figure out the 
other graphlab types easily.
\code
typedef graphlab::types<graph_type> gl;
\endcode
Rather than needing to directly instantiate template interfaces like:
\code
graphlab::iscope< graphlab::<vertex_data, edge_data> >* scope;
\endcode
we can use the simpler syntax:
\code
gl::iscope* scope;
\endcode

\section detailed_example_graph_datastructure The Graph Data Structure

The next step in constructing a GraphLab program is to construct the actual graph. 
The core graph datastructure is documented here: graphlab::graph .
To simplify the presentation we will define a function which takes a reference to a 
graph and populates the graph. Later we will show how to construct the empty graph. 
The init_graph function takes an additional argument which describes the number of 
coins along each dimension of the grid.
\code 
void init_graph(graph_type& g, size_t dim) {
\endcode

We first create create d * d vertices. We use the graph's 
\ref graphlab::graph::add_vertex "gl::vertex_id_t add_vertex(vertex_data)"
method which takes the vertex data as input and returns the vertex 
id of the new vertex. The ids are guaranteed to be sequentially numbered. 
The graph data structure behaves like an STL container and stores the vertex data by value.

\code
  for (size_t i = 0; i < dim * dim; ++i) { 
    // create the vertex data, randomizing the color 
    vertex_data vdata; 
    vdata.numflips = 0; 
    // Flip a uniform coin to obtain the initial color 
    if (gl::random::rand_int(1) == 1) vdata.color = true; 
    else vdata.color = false; 
    // create the vertex 
    g.add_vertex(vdata); 
  }
\endcode

Now we add all the edges. To add edges we use the graph's 
\ref graphlab::graph::add_edge "gl::edge_id_t add_edge(src, target, edgedata)"
function which creates an edge from src to target with the
edge data given by edgedata. The add_edge function then returns the 
id of the new edge. The ids are guaranteed to be sequentially numbered. 
GraphLab does NOT support duplicated edges, and currently has no facilities 
for checking for accidental duplicated edge insertions at the graph construction
stage. (It is quite costly to do so) Any duplicated edges will result in an 
assertion failure at the later finalize stage. Furthermore, the current version 
does not support vertex or edge removal. These constraints are imposed to enable 
the efficient construction of massive graphs while retaining fast look-up. For 
more details about the graph data-structure see \ref graphlab::graph
\code
edge_data edata; 
for (size_t i = 0; i < dim; ++i) { 
  for (size_t j = 0; j < dim - 1; ++j) { 
    // add the horizontal edges in both directions 
    g.add_edge(dim * i + j, dim * i + j + 1, edata); 
    g.add_edge(dim * i + j + 1, dim * i + j, edata); 
    // add the vertical edges in both directions 
    g.add_edge(dim * j + i, dim * (j + 1) + i, edata); 
    g.add_edge(dim * (j + 1) + i, dim * j + i, edata); 
  } 
}
\endcode
The above block of code connects all the vertices in the 
grid pattern illustrated in Figure 1(a). Now that the graph is fully 
constructed, we need to call \ref graphlab::graph::finalize graph.finalize().

\code
  g.finalize(); 
} // end of init_graph function
\endcode

The finalize function reorders the vertex adjacency tables so that 
the \ref graphlab::graph::in_edge_ids "in_edge_ids(vertex_id_t)" returns edges in order of the source 
vertex id and the \ref graphlab::graph::out_edge_ids "out_edge_ids(vertex_id_t)" returns edges in order of 
the target vertex id. This sorting also enables O(log(degree)) edge retrieval. 





\section detailed_example_update_function The Update Function

Now we define the update function which represents the basic block of 
computation in this program. The update function is applied to each vertex 
and has read/write access to the data at that vertex, as well as all adjacent 
edges and vertices. You may specify more than one update function, but we only 
need one for this application. Lets first outline the computation that will 
occur in the update function. Below is a list of the steps:

\li \b 1. Compute the sums of the red and black neighbors.
\li \b 2. Draw the new assignment from Bernoulli(prop_red)
\li \b 3. Count a change if a new assignment was drawn.
\li \b 4. Decide whether to reschedule this vertex and any of its neighbors: \n
    \b a) The color was not changed and so neighbors do not need to be re-updated. Otherwise update neighbors. \n
    \b b) If this was a deterministic decision (all neighbors had the same color) we do not need to reschedule this vertex. Otherwise update this vertex.

Each update function call has the option of inserting new 
tasks into the scheduler: in this case, its self and its neighbors. 
The algorithm terminates when there are no tasks remaining. There are 
other methods for terminating execution, such as registering a termination 
valuator with the engine, but we are not going to describe that here. We now 
present the entire update function and then discuss each part individually.

\code
void update_function(gl::iscope& scope,
                     gl::icallback& scheduler) {
  // Get a reference to the vertex data an in edges
  vertex_data& curvdata = scope.vertex_data();
  graphlab::edge_list in_edges = scope.in_edge_ids();
  // Count the number of red neighbors
  size_t num_red_neighbors = 0;  
  for (size_t i = 0; i < in_edges.size(); ++i) {
    // eid is the current edge id
    size_t eid = in_edges[i];    
    size_t sourcev = scope.source(eid);
    const vertex_data& nbrvertex = scope.neighbor_vertex_data(sourcev);
    if (nbrvertex.color) ++num_red_neighbors;
  }
  // get the total number of neighbors we have
  size_t num_neighbors = in_edges.size();
  // Draw the new color
  bool new_color =
    gl::random::rand01() < (double(num_red_neighbors) / num_neighbors);
  // Determine if the draw was deterministic
  bool is_deterministic =
    num_neighbors == num_red_neighbors || num_red_neighbors == 0;
  // see if I flip and update the current vertex data.
  bool color_changed = new_color != curvdata.color;
  if (color_changed) ++curvdata.numflips;
  // Assign the new color 
  curvdata.color = new_color;
  // If I flipped, all my neighbors could be affected, loop through 
  // all my neighboring vertices and add them as tasks. 
  if (color_changed) {
    for (size_t i = 0; i < in_edges.size(); ++i) {
      size_t sourcev = scope.source(in_edges[i]);
      scheduler.add_task(gl::update_task(sourcev, update_function),
                        1.0);
    }
  }
  // Reschedule myself if this was not a deterministic draw
  if (is_deterministic == false) {
    scheduler.add_task(gl::update_task(scope.vertex(), update_function),
                        1.0);
  }
}
\endcode

All update functions must have the following form:

\code
void update_function(gl::iscope& scope, 
                     gl::icallback& scheduler) {
\endcode

The parameters are described here:

\par Scope:
    The scope provides access to a local neighborhood of a graph and has the type \ref graphlab::iscope .
    The scope is centered on a particular vertex, scope.vertex(), and includes 
    all adjacent edges and vertices. All vertices are identified by an unsigned 
    integer type vertex_id_t, and all edges are similarly identified by an 
    unsigned integer type edge_id_t. GraphLab guarantees that all vertices are 
    sequentially numbered from 0 (so the largest vertex id is <tt>|num_vertices| - 1</tt>), 
    and similarly for edges. All edges are directed. The scope provides methods 
    to access the data associated with the vertex its neighbors and any inbound or outbound edges. 
    
\par Scheduler:
    There are three basic types of schedulers. The first type consists only of the 
    synchronous scheduler which is automatically used with the synchronous engine. 
    The second type of schedulers are the static asynchronous scheduler like the 
    chromatic scheduler. These static schedulers execute a fixed static schedule 
    until some termination condition is reached. None of the first two types of 
    schedulers used the scheduler callback.
    The last class of schedules are the task schedulers, which enable dynamic 
    computation and can receive new tasks (update function, vertex pairs) from 
    the update functions by using the icallback interface. The schedulers include 
    fifo, multiqueue_fifo, priority, multiqueue_priority and clustered_priority. 

First we get a mutable reference to the vertex data on this vertex. 

\code 
  vertex_data& curvdata = scope.vertex_data();
\endcode

Note that the function \ref graphlab::iscope::vertex_data "scope.vertex_data()" returns a 
reference to the vertex data. Modifications to this reference will directly modify the data 
on the graph. Then we get a constant reference to the vector of edge ids for all edges inbound 
to this vertex: 

\code
  const std::vector<gl::edge_id_t>& in_edges = scope.in_edge_ids();
\endcode

We now compute the total number of red neighbors. This is done by looping through all my 
neighboring vertices and counting the number of red vertices. To do this I have to look 
at my neighboring vertices data. The \ref graphlab::iscope::neighbor_vertex_data "neighbor_vertex_data(vid)" 
function allow me to read  the vertex data of a vertex adjacent to the current vertex. If we have edge data, the 
function \ref graphlab::iscope::edge_data "edge_data(eid)" will return a reference to the edge 
data on the edge eid. Since  I am not going to need to change this data, I can just grab a const reference. 
You should  always try to use const references whenever you know that you will definitely not be 
changing the data, since GraphLab could make additional optimizations for "read-only" 
operations. Similarly, you should never cast a constant reference to a regular reference. 
Modifications to constant references have undefined behavior.

\code
  // Count the number of red neighbors
  size_t num_red_neighbors = 0;  
  for (size_t i = 0; i < in_edges.size(); ++i) {
    size_t eid = in_edges[i];    
    size_t sourcev = scope.source(eid);
    const vertex_data& nbrvertex = scope.neighbor_vertex_data(sourcev);
    if (nbrvertex.color) ++num_red_neighbors;
  }
\endcode

Now we can decide on the new color of the vertex. We either match the 
majority color, or in the case of no majority color, we flip a coin. 
Extensive random number support is provided through gl::random. 
random::rand01() provides a random floating point number 
between 0 and 1. There are a number of other distribution based generators available in \ref random

\code
  // Draw the new color
  bool new_color =
    gl::random::rand01() < (double(num_red_neighbors) / num_neighbors);
\endcode

Now we have decided on the new color of the vertex, we can go ahead and update the 
value of the vertex. Once again recall that curvdata is a reference to the actual
graph data. Therefore we can just modify it directly. Here we will track whether the 
flip was deterministic and whether a new color was chosen.

\code
  // Determine if the draw was deterministic
  bool is_deterministic =
    num_neighbors == num_red_neighbors || num_red_neighbors == 0;
  // see if I flip and update the current vertex data.
  bool color_changed = new_color != curvdata.color;
  if (color_changed) ++curvdata.numflips;
  // Assign the new color 
  curvdata.color = new_color;
\endcode

Now for the task creation algorithm. There are 2 basic cases:

\li <b>Case 1:</b> \n
    If I did not flip, then my neighbors will not be affected, it will be as 
    if this vertex was never updated at all. We therefore do not need to update my neighbors. 
\li <b>Case 2:</b> \n
    If I flipped, all my neighbors could be affected, therefore loop through 
    all my neighboring vertices and add them as tasks. To add a task, I call 
    \ref graphlab::icallback::add_task "scheduler.add_task(gl::update_task, priority)". 
    The gl::update_task object takes a vertex id, and the update function to execute on. The priority 
    argument is the priority of the task. This is only used by the priority schedulers. 
    This value should be strictly > 0. In this demo app, we don't really care about the 
    priority, so we will just set it to 1.0. 
  
\code
  if (color_changed) {
    for (size_t i = 0; i < in_edges.size(); ++i) {
      size_t sourcev = scope.source(in_edges[i]);
      scheduler.add_task(gl::update_task(sourcev, update_function),
                        1.0);
    }
  }
\endcode
Now, there is another special case. If flipped myself on a random number, 
then I could switch colors when updating myself again. Therefore I should 
try again and update myself again in the future 
\code
  // Reschedule myself if this was not a deterministic draw
  if (is_deterministic == false) {
    scheduler.add_task(gl::update_task(scope.vertex(), update_function),
                        1.0);
  }
} // end of update_function
\endcode

\section detailed_example_shared_variables Shared Variables
The shared variable system serves two roles. First, it provides 
access to data which is globally accessible through all update functions. 
Second, the shared variables provides the capability to compute aggregations 
of all graph data. The first capability may not appear to be useful in the 
shared memory setting since one could simply define global variables. 
However, using the shared data manager, allows GraphLab to manage these
"global variables" in a platform independent fashion; such as in the 
distributed setting and will lead to more portable code. Additionally, the Shared variables
system use an RCU mechanism to ensure safe access to shared data.

For this application, we are interested in having an incremental 
counter which provides the total number of flips executed so far, 
as well as computing the proportion of red vertices in the graph. 
We can achieve these tasks using the shared data Sync mechanism.

The Sync mechanism allows you to build a 'Fold / Reduce' operation 
across all the vertices in the graph, and store the results in a shared variable.

\li the total number of vertices (constant)
\li red vertex proportion (synced)
\li the total number of flips (synced)

We will therefore define the following 3 shared variables: 
\code
gl::glshared_const<size_t> NUM_VERTICES;
gl::glshared<double> RED_PROPORTION;
gl::glshared<size_t> NUM_FLIPS;
\endcode

\subsection detailed_example_num_vertices Num Vertices

The number of vertices in the graph is a constant, 
and will not be changed through out execution of the GraphLab program. 
We simply just set its value using the set() function.
\code
NUM_VERTICES.set(DIM * DIM);
\endcode
To get the value of a shared variable:
\code
size_t numvertices = NUM_VERTICES.get_val();
\endcode


\subsection detailed_example_red_proportion Red Proportion

A sync is defined by a pair of functions, a reducer, and an apply The reducer is 
exactly a fold over all the vertices, and the apply takes the final value at the 
end of the reduce, and performs whatever transformation it needs, before writing 
it into the Shared variable. For instance, an L2 sum can be computed by having the 
reducer add the squares of the values at each vertex, then the apply function 
performs the square root.

We will use this to implement the RED_PROPORTION sync. The way we will 
implement this is to use the reducer to count the number of Red vertices. 
The apply function will then divide the result by the value in the NUM_VERTICES 
table entry. The reducer is a function is of the form:
\code
void reduce_red_proportion( gl::iscope& scope, 
                            graphlab::any& accumulator) {
\endcode

\par scope
    The scope on the vertex we are currently accessing. 
\par accumulator
    The input and output of the fold/reduce operation.

In this reducer, we will simply increment the accumulator if the color of the vertex is red.
\code
  if (scope.vertex_data().color) accumulator.as<double>()++; 
}
\endcode

The apply function takes the followng 2 parameters

\code
void apply_red_proportion(graphlab::any& current_data, 
                          const graphlab::any& new_data) {
\endcode

\par current_data
    The current (old) value of the shared variable. Overwriting this will 
    update the value of the shared variable. The type of the data matches 
    the type of the shared variable.
\par new_data
    The result of the reduce operation.

The reduced result in new_data, will be the count of the number of red vertices. 
We can get the total number of vertices from NUM_VERTICES, and compute the 
proportion of red vertices by dividing the number of red vertices by the total 
number of vertices.
\code
  size_t numvertices = NUM_VERTICES.get(); 
  // new_data is the reduced result, which is the number 
  // of red vertices 
  double numred = new_data.as<double>(); 
  // compute the proportion 
  double proportion = numred / numvertices; 
  // here we can output something as a progress monitor 
  std::cout << "Red Proportion: " << proportion << std::endl; 
  // write the final result into the shared variable 
  current_data.as<double> = (double)proportion; 
}
\endcode

Now that both reduce and apply functions have been defined, we can create the
sync, by calling in the beginning of the program
\code
core.set_sync(RED_PROPORTION, 
              reduce_red_proportion, 
              apply_red_proportion, 
              double(0), 
              100);
\endcode
Here we have set the reduction proportion key to store the 
result of running reduce_red_proportion to fold over all the 
vertices with starting value double(0) and result applied by 
using apply_red_proprotion. The operation will be applied 
approximately every 100 updates. 


\subsection detailed_example_num_flips Num Flips

GraphLab provides a number of predefined syncing operations which allow simple 
reductions / applies to be implemented very quickly. For instance, computing 
sums, sums of squares, etc. We will implement the NUM_FLIPS entry using one 
of these predefined operations. Since the vertex data could be any arbitrary type, 
the predefined operations typically require the user to provide a simple function 
which extracts the information of interest from the vertex data. In this case, 
we are interested the numflips field.
\code
size_t get_flip(const vertex_data& v) { 
  return v.numflips; 
}
\endcode

To create the sync, we use the set_sync function as well, 
but using functions from \ref graphlab::glshared_sync_ops and \ref graphlab::glshared_apply_ops.
In this case, our reduction function is a simply "sum", while our apply 
function should do nothing more than copy the result of the reduction 
into the shared variable.
\code
core.set_sync(NUM_FLIPS, 
              gl::glshared_sync_ops::sum<size_t, get_flip>, 
              gl::glshared_apply_ops::identity<size_t>, 
              size_t(0), 
              100);
\endcode

\subsection detailed_example_sync_together Putting It Together
We can now write a simple init_shared_data function to take 
a reference to a core and initialize all the shared variables.

\code
void init_shared_data(gl::core& core, size_t dim) { 
  NUM_VERTICES.set(dim * dim); 
  core.set_sync(RED_PROPORTION, 
                reduce_red_proportion, 
                apply_red_proportion, 
                double(0), 
                100); 
  core.set_sync(NUM_FLIPS, 
                gl::glshared_sync_ops::sum<size_t, get_flip>, 
                gl::glshared_apply_ops::identity<size_t>, 
                size_t(0), 
                100); 
} 
\endcode

\section detailed_example_merges Merge Function

A merge operation can also be provided which will allow the 
sync operation to be parallelized. This is highly recommended.

If a merge function is defined for a sync, the set of vertices will be 
partitioned into a collection of disjoint sets. The sync function is then 
performed on each set in parallel, producing a number of partial results. 
These partial results are then combined with the merge function.
Finally, the apply function is executed on the final result and written to the shared variable.

\subsection detailed_example_merge_red Merging Red Proportion
Since the intermediate results of the sync is simply a count of the number of red vertices
seen so far, the merge is then just a sum.
\code
void merge_red_proportion(graphlab::any& target, 
                          const graphlab::any& source) {
  target.as<double>() += source.as<double>();
}
\endcode


To sync using the merge operation, the set_sync is called using:
\code
core.set_sync(RED_PROPORTION, 
              reduce_red_proportion, 
              apply_red_proportion, 
              double(0), 
              100,
              merge_red_proportion);
\endcode

\subsection detailed_example_merge_red Merging Num Flips
Similarly, a number pre-defined merge functions are defined in \ref graphlab::glshared_merge_ops ,
and the NUM_FLIPS sync can be defined in the same way.
\code
core.set_sync(NUM_FLIPS, 
              gl::glshared_sync_ops::sum<size_t, get_flip>, 
              gl::glshared_apply_ops::identity<size_t>, 
              size_t(0), 
              100,
              gl::glshared_merge_ops::sum<size_t>);
\endcode


\section detailed_example_main The Main 

The Main function where everything begins. Here, we will 
demonstrate the minimal code needed to start a GraphLab job 
using all the parts we defined above. We will use the GraphLab 
command line tools to setup a basic engine using command line options. 
We first present the complete main and then discuss each of the parts.

\code
int main(int argc,  char *argv[]) {
  // Parse command line options
  graphlab::command_line_options opts;
  size_t dimensions = 20;
  opts.attach_option("dim",
                     &dimensions, dimensions,
                     "the dimension of the grid");
  opts.scheduler_type = "fifo";
  opts.scope_type = "edge";
  if(!opts.parse(argc, argv)) return EXIT_FAILURE;
  // Create the core which contains the graph and engine
  gl::core glcore;
  // Initialize engine with command line options
  glcore.set_engine_options(opts);
    // Initialize the the data structures  
  init_graph(glcore.graph(), dimensions);
  init_shared_data(glcore, dimensions);
  // Add all starting tasks
  glcore.add_task_to_all(update_function, 1.0);  
  // Run the graphlab engine 
  double runtime = glcore.start();
  
  // Output the results
  std::cout << "Completed in " << runtime << " seconds" << std::endl;
  glcore.sync_now(NUM_FLIPS);
  glcore.sync_now(RED_PROPORTION);
  // now we can look the values using the get() function
  size_t numberofflips = NUM_FLIPS.get_val();
  double redprop = RED_PROPORTION.get_val();
  std::cout << "Number of flips: " <<  numberofflips << std::endl;
  std::cout << "Red prop: " << redprop << std::endl;
  // output the graph
  size_t ctr = 0;
  for (size_t i = 0;i < dimensions; ++i) {
    for (size_t j = 0;j < dimensions; ++j) {
      std::cout << size_t(glcore.graph().vertex_data(ctr).color) << " ";
      ++ctr;
    }
    std::cout << std::endl;
  }
}
\endcode

Since the GraphLab engine can take many options we have built-in some 
command line parsing tools. In this program we add an additional command 
line argument "dim" which specifies the size of the grid. We set the default 
value to 20. In addition we set the default scheduler to be "fifo" and the 
default scope type to be "edge". If the user provides different values than 
the defaults will be replaced. The call to \ref graphlab::command_line_options::parse "opts.parse(argc,argv)" invokes the p
arser and fills in the fields.

The \ref graphlab::core "gl::core" object bundles an empty graph, and engine configuration into 
a single object. The \ref graphlab::core::set_engine_options "core.set_engine_options(opts)" takes the engine 
options from the command line and uses them to configure the internal engine. 
The graph can be retrieved from the glcore by calling the \ref graphlab::core::graph "glcore.graph()" function. 
The \ref graphlab::core::add_task_to_all "glcore.add_task_to_all" function adds an update task to each vertex with the 
desired update function and priority value.

Finally, the engine is run by calling \ref graphlab::core::start "glcore.start()" which runs until their
are no more tasks remaining or until a termination condition is reached.

\section detailed_example_comand_line Running On the Command Line
After compilation, 
\verbatim
./demo --help
\endverbatim
will produce a list of all the available options.
\verbatim
GraphLab program.:
  --help                  Print this help message.
  --dim arg (=20)         the dimension of the grid
  --ncpus arg (=2)        Number of cpus to use.
  --engine arg (=async)   Options are {async, async_sim, synchronous}
  --affinities arg (=0)   Enable forced assignment of threads to cpus
  --schedyield arg (=1)   Enable yielding when threads conflict in the 
                          scheduler.
  --scope arg (=edge)     Options are {none, vertex, edge, full}
  --metrics arg (=basic)  Options are {none, basic, file, html}
  --schedhelp arg         Display help for a particular scheduler.
  --scheduler arg (=fifo) Supported schedulers are: chromatic, sweep, fifo, 
                          priority, multiqueue_fifo, multiqueue_priority, 
                          splash, round_robin, clustered_priority, sampling. 
                          Too see options for each scheduler, run the program 
                          with the option ---schedhelp=[scheduler_name]
\endverbatim
Observe that the <tt>dim</tt> option defined by the \ref graphlab::command_line_options::attach_option "attach_option()"
call in the main function appears as an available option.

The GraphLab command line permit quite flexible manipulation of the scheduler capabilities
and options through the command line using the <tt>--scheduler arg</tt> options.
Running <tt>--schedhelp</tt> displays all the available scheduler options.

When this tutorial was written, the output of running <tt> ./demo --schedhelp </tt> is

\verbatim
chromatic scheduler
--------------------------------------------------
a scheduler which performs #iterations sweeps of
the graph using a graph color ordering.

Options: 
max_iterations = [integer, default = 0]
update_function = [update_function_type,default = set on add_task]

sweep scheduler
--------------------------------------------------
very fast dynamic scheduler. Scans all vertices in
sequence, running all update tasks on each vertex
evaluated.

Options: 
ordering = [string: linear/permute, default=linear]

fifo scheduler
--------------------------------------------------
Standard FIFO task queue, poor parallelism, but
task evaluation sequence is highly predictable.
Useful for debugging and testing.

Options: 

priority scheduler
--------------------------------------------------
Standard Priority queue, poor parallelism, but
task evaluation sequence is highly predictable.
Useful for debugging

Options: 

multiqueue_fifo scheduler
--------------------------------------------------
One or more FIFO task queues is assigned to each
processor, where the queues are stochastically
load balanced. Like the fifo scheduler, but less
predictable, and much faster.

Options: 

multiqueue_priority scheduler
--------------------------------------------------
One or more Priority task queues is assigned to
each processor, where the queues are
stochastically load balanced. Like the priority
scheduler, but less predictable, and much faster.

Options: 

splash scheduler
--------------------------------------------------
Similar to the priority queue scheduler, but
allows for only one update function. Updates are
evaluted in a "splash" ordering

Options: 
splash_size = [integer, default = 100]
update_function = [update_function_type,default = set on add_task_to_all]

round_robin scheduler
--------------------------------------------------
Loops over a sequence of tasks repeatedly for #
iterations.

Options: 
max_iterations = [integer, default = 0]
start_vertex = [integer, default = 0]

clustered_priority scheduler
--------------------------------------------------
Like the priority scheduler, but groups vertices
into clusters where the entire cluster has a
single priority

Options: 
partition_method = [string: metis/random/bfs, default=metis]
vertices_per_partition = [integer, default = 100]

sampling scheduler
--------------------------------------------------
A scheduler which samples vertices to update based
on a multinomial probability which can be updated
dynamically.

\endverbatim

This lists all the available schedulers and the available options for each scheduler.
For instance: to run the demo process using the sweep scheduler with a randomly permuted ordering:

\verbatim
./demo --scheduler="sweep(ordering=permute)"
\endverbatim

Additional scheduler options are seperated with a comma. For instance to run with the clustered
priority scheduler using random partitioning and 50 vertices per partition:

\verbatim
./demo --scheduler="clustered_priority(partition_method=random,vertices_per_partition=50)"
\endverbatim


\section detailed_example_engine_options Comments on Engine Options
\subsection detailed_example_scope_model Scope Model
When an update function is executed on a vertex, 
it can access all graph data on adjacent edges and adjacent 
vertices. The different scoping consistency models provide 
different data consistency guarantees when accessing graph data. 
There are three scoping models, vertex, edge, and full.

\par Vertex Consistency
Vertex consistency is the weakest consistency model, and also the fastest 
(lowest contention). The vertex consistency model only guarantees that the 
Update Function can read and write to the current vertex without experiencing 
data races. Reading or writing to adjacent edges or adjacent vertices could 
result in inconsistent data. Data on adjacent edges and vertices may also 
change between consecutive reads within a single Update Function call.
\par Edge Consistency
The edge consistency model guarantees that the Update Function can read and 
write to the current vertex as all as all adjacent edges without experiencing 
data races. In addition, the Update Function can also made consistent reads 
from adjacent vertices.
\par Full Consistency
The full consistency model guarantees that the Update Function can read and 
write to the current vertex, as well as all adjacent edges and vertices in 
a consistent fashion. This model experiences the highest amount of contention 
and provides the lowest level of parallelism

\subsection detailed_example_consistency Choosing a consistency model

The user should try to pick the lowest consistency model which satisfies 
the needs of the algorithm. For instance, in this demo application, since 
the update function only requires reading of neighboring vertex data, 
the edge_consistency model is guaranteed to have sequential consistency, 
and the algorithm is therefore guaranteed to be correct (assuming GraphLab is bug-free) 
if executed with the edge consistency model or the full consistency model.

Note that the sync operation is guaranteed to be sequentially consistent

\subsection detailed_example_scheduler Scheduler Type

GraphLab provides eight schedulers. The Synchronous scheduler, 
the Round-robin scheduler, five task schedulers, and the Splash scheduler.

All four task schedulers behave similarly as in the demo application, 
but each have different set of scheduling guarantees.

\par FIFO scheduler (fifo): 
  Implements a strict single FIFO queue of tasks
\par Multiqueue FIFO scheduler (multiqueue_fifo): 
  Uses multiple load balanced FIFO queues to decrease contention. This tends to have better performance over FIFO, but loses the "First-in-first-out" guarantees.
\par Sweep scheduler (sweep): 
  partitions the vertices among the processors. Each processor than loops through all the vertices in its partition, executing all tasks encountered.
\par Priority scheduler (priority): 
  Implements a priority queue over tasks. Executes tasks in priority order.
\par Multiqueue Priority scheduler (multiqueue_priority): 
  Like Multiqueue FIFO scheduler but with prioritized tasks. Weaker priority guarantees but better performance.
\par Clustered Priority scheduler (clustered_priority): 
  partitions the graph into a collection of subgraphs, and builds a priority queue over subgraphs. tasks within a subgraph are executed in arbitrary order. The partitioning methods are "metis","bfs" and "random". Metis provides the best partitioning, but could be extremely costly for large graphs.

*/

