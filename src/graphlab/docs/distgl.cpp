/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


/**

\page distributed_gl Using Distributed GraphLab
Distributed GraphLab is necessarily more complex than the shared memory
version, though much effort has gone into making it as useable as possible.

The distributed GraphLab implementation is built on top of the \ref RPC "RPC"
framework, and at least a cursory glance at \ref Spawning
as well \ref graphlab::distributed_control might be useful.

\li \ref page_distributed_graph_creation
\li \ref page_distributed_graphlab


\page page_distributed_graph_creation Distributed Graph Creation
The goal of distributed GraphLab is to both handle graphs that exceed the memory
capacity, and to make use of more processors than what is available on a 
single node. Since the number of machines available may 
vary with performance demands, the graph representation must 
support loading of the data graph on varying sized cluster deployments. 
To resolve these issues, we developed an on-disk graph representation based 
on a two-phased partitioning scheme. 

\section sec_disk_graph_atom Atom partitioning
The initial user graph must first be over-partitioned into 
more parts than the largest cluster deployment you will ever use for the data.
Each of these parts is called an <b> atom </b>. An atom not just contains the
vertices and edges within its partition, but also contains information about the
vertices and edges immediately adjacent to the partition (also called the 
<b> ghost </b> of the partition.) 

A collection of atoms which make up a complete graph is called a <b> disk graph </b>.
The atom file names are of the form "[basename].0", "[basename].1",
"[basename].2", etc. Each atom is stored a binary format 
Kyoto Cabinet (http://fallabs.com/kyotocabinet/) hash table. 
The user should almost never have to interact with the atoms directly, but the format
is documented briefly in graphlab::disk_atom. 

A proper disk graph must also include an <b> atom index file </b> with 
filename "[basename].idx". The atom index is an additional human readable/modifiable 
text file which provides some basic
information about the graph (number of vertices, number of edges), the adjacency 
structure of the atoms, and most importantly, the <b> file location / URL </b>
of the atom files themselves. (Currently only the file:// URL protocol is supported,
but future enhancements might support additional protocols.)

\note In cases where the atom index is missing, the \ref graphlab::disk_graph does
provide a constructor to load a disk graph without an atom index. After which 
finalization will re-generate the atom index.


\section sec_disk_graph_creation Graph Creation

There are three different methods of creating a disk graph, targeted at three different 
use cases.

\li \ref sec_memory_to_disk_graph Graph is small enough to <b>fit in memory</b> on one machine.
This case is common when the user is testing a distributed GraphLab implementation, 
or when computation time is large compared to the graph size.

\li \ref sec_disk_graph Graph is small enough to be <b>created on one machine</b>.
This case is the middle ground between the methods above and below this. The graph is small
enough that it can be written to disk in reasonable time on one machine.

\li \ref sec_distributed_disk_graph Graph is too large to create on one machine. Distributed
processing is necessary to even construct the graph in reasonable time. This is the most complex
case and requires more user effort than the above two.

\subsection sec_memory_to_disk_graph Conversion of In Memory Graph to Disk Graph
This is the simplest case and the user only has to instantiate a regular in-memory graph
using graphlab::graph. Following which an atom partitioning could be generated
using any one of the partition methods in graphlab::graph. 
Calling graphlab::graph_partition_to_atomindex will convert the graph to a disk graph.

\subsection sec_disk_graph Direct Creation of Disk Graph
\copydoc graphlab::disk_graph

\subsection sec_distributed_disk_graph Distributed Graph Creation
A MapReduce style distributed graph construction method is provided
by graphlab::mr_disk_graph_construction. The basic requirement is that all the machines
have a shared file storage (NFS / SMBFS / etc).

The user first has to create a subclass of the \ref graphlab::igraph_constructor
which simply provides a iteration mechanism over an arbitrary subset of the vertices/edges
of the graph. Some restrictions on the graph constructor are necessary since it must be
"distributed-aware".

After which \ref graphlab::mr_disk_graph_construction will use the graph constructor to
perform a highly parallel/distributed construction of the disk graph. 

This is best illustrated using an example
\ref distributed_dg_construction_test.cpp



*/

/**
\page page_distributed_graphlab Distributed GraphLab
To work with distributed GraphLab, MPI is needed. We have tested with MPICH2 and the MPI
distribution that comes with OS X. However, other distributions of MPI (Open MPI) should
work as well.

Distributed GraphLab is functionally quite similar to the regular GraphLab. But due to the
law of leaky abstractions (http://www.joelonsoftware.com/articles/LeakyAbstractions.html), 
there are some issues the user will have to take into consideration.

Firstly, the graph cannot be created on the fly, and must be created before-hand and loaded from 
an atom index file ( \ref page_distributed_graph_creation ).

Next, since multiple instances of the program is started across the network, some amount of 
coordination is needed to ensure that all instances behave identically. We strongly encourage
the use of graphlab::core (which is the distributed setting is graphlab::distributed_core) since
that simplifies the amount of management needed. 

\attention The general rule is that any operation which affect operation on a global scale
(for instance, a configuration setting) should be called by all instances simultaneously.

For users who have read the shared memory detailed_example, or are already familiar with 
the shared memory GraphLab, can simply read on. 
Otherwise here is a distributed version of the shared memory detailed_example: 
\ref distributed_detailed_example

Depending on the complexity of your program, moving from the
shared memory version to the distributed version can be quite mechanical.
The key differences are listed here.

An easy way to see what changes you have to make to move a shared memory version
to a distributed version is to take a diff between tests/demo/demo.cpp and tests/dist_demo/dist_demo.cpp

\section sec_distributed_types Types
Instead of including <tt>graphlab.hpp</tt>, you should include <tt>distributed_graphlab.hpp</tt>
Similarly, the distributed version core is called <tt>dgraphlab::istributed_core</tt> and the distributed
version of the types struct is <tt>graphlab::distributed_types</tt>.

We recommend using the <tt>graphlabb::distributed_types</tt> system,
\code
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::distributed_types<graph_type> gl;
\endcode
since it this will allow you to "copy and paste" code from your shared memory implementation
easily while the type system changes the back-end implementations.

For instance, your update function should not use <tt>graphlab::edge_list</tt>, but should
use <tt>gl::edge_list</tt> since the exact type is different for the shared memory graph
and for the distributed graph.

Similarly, <tt>gl::distributed_core</tt> should be used instead of <tt>core</tt>  and 
<Tt>gl::distributed_glshared</tt> should be used instead of <tt>glshared</tt>

\section sec_distributed_use_disk_graph Use The Disk Graph
Unlike the shared memory version, the graph cannot be created on the fly. The user must first
create a disk graph, then construct the distributed_core / distributed_graph from the disk graph.

\section sec_distributed_build_parallel_execution Distributed Execution
The user must constantly keep in mind that in the distributed setting, there are a collection
of "independent" copies of the program being executed, all communicating through MPI or
through the \ref RPC "GraphLab RPC" system. 

All distributed_core / engine functions generally require all machines to execute the exact
same code path, running the same functions at the same time. For instance, the following code:
\code
  core.set_sync(...);
  core.set_engine_options(clopts);
  core.build_engine();
\endcode
should be executed by all machines in the same order at the same time. Some of these functions
(for instance <tt>core.build_engine()</tt> ) has an internal distributed barrier which ensures that
all machines must reach the <tt>build_engine()</tt> call before execution can proceed.

There are a few exceptions to this rule (such as <tt>add_task()</tt>) and these are documented.

\section sec_distributed_graph_limitations Distributed Graph Limitations
The regular graph has functions which permit obtaining vertex and edge data by reference.
The distributed_graph has similar functions but for obvious reasons, it is restricted to 
vertex and edge data "owned" by the local partition. The alternate functions
get_vertex_data() / set_vertex_data() as well as the edge varieties should be used instead.
These return by value and can get/set vertex and edge data from across the network. 

The user should keep in mind that these remote network gets/sets are inherently slow
and is where the Law of Leaky Abstractions slip through. If a large number of vertex/edge reads
and writes are needed, it is recommended that the user make use of the \ref RPC "RPC" operations
to gather/scatter the information in bulk instead of accessing the data one at a time. Alternatively,
the graph can be saved to disk and loaded using one machine making use of the disk_graph 
class.

\section sec_distributed_graphlab_build_the_engine Build The Engine
The shared memory core automatically destroys and reconstructs the engine as options change.
This is unfortunately difficult to coordinate in the distributed setting. Therefore the user
must explicitly construct the engine once all options are set through the \ref graphlab::distributed_core::build_engine()
function. 

Once the engine has been constructed, it can be used repeatedly. However, engine options
can no longer be modified.

\section sec_distributed_graphlab_command_line_options Distributed Command Line Options
The graphlab::command_line_options object automatically injects options for the shared
memory engines. To get options for the distributed engine, 
\code
  opts.use_distributed_options();
\endcode
must be called prior to calling <tt>opts.parse()</tt>

\section sec_distributed_graphlab_engine_types Distributed Engine Types
The user will notice that the new command line options include a <tt>--engine</tt> option.
Distributed Graphlab implements two engines with very different characteristics. One is
called the chromatic engine (<tt>--engine=dist_chromatic</tt>), the other is called the locking engine
<tt>--engine=dist_locking</tt>. The chromatic engine is the default.

\subsection sec_distributed_graphlab_chromatic Chromatic Engine
The chromatic engine operates by associating an update_task with each vertex. It then makes
use of the graph coloring to run update tasks in parallel. For instance if my graph has 3 colors,
it will run all update tasks on vertices of color 0 in parallel. Stop, then run all tasks
on vertices of color 1 in parallel, etc.

The chromatic engine therefore does not use a scheduler since it essentially has a built in 
chromatic scheduler. However, it assumes that <b>The graph coloring is valid for the 
scope type requested</b>. In other words, if an edge scope is requested, then the graph coloring
must ensure that neighboring vertices do not have the same color. If a full scope is requested,
then the graph coloring must ensure that vertices one hop away from each other do not have the 
same color. If coloring constraints are violated, the engine will still run properly, though
consistency guarantees do not hold.

\subsubsection chromatic_engine_options Engine Options
The max iteration count and whether to permute 
the vertex update order can be set. There are two engine options, "max_iterations=N"
and "randomize_schedule=0 or 1" and these can be set on the command line:
\verbatim
--engine="dist_chromatic(max_iterations=10,randomize_schedule=1)"
\endverbatim

\subsubsection chromatic_engine_limitations Limitations
The chromatic engine is therefore somewhat limited in its scheduling capabilities
(only supporting a chromatic schedule), supports only one update function type, and requires
a valid graph coloring. However, it benefits from having a much smaller overhead than the
other more general Locking engine. It is therefore the default engine.

\subsection sec_distributed_graphlab_locking Locking Engine
The locking engine is the generalization of the shared memory implementation to the distributed
setting. It makes use of distributed locking to ensure consistency, and uses a variety of pipelining
methods to keep a large number of operations in flight at any time. It can work with a subset of the 
shared memory schedulers.

The locking engine is much more general and operates almost identically to the shared memory
engines. However, the implementation complexity leads to some loss of efficiency.

\subsubsection locking_engine_options Engine Options
The locking engine takes one option: max_deferred_tasks_per_node which is the maximum number
of update tasks to keep in flight in the locking pipeline. The default value is 1000.
This can be set on the command line using:
\verbatim
--engine=dist_locking(max_deferred_tasks_per_node=10)"
\endverbatim

\section sec_distributed_graphlab_not_yet_implemented Incomplete Implementations
<tt>distributed_glshared_const</tt> is not yet implemented. The user can "emulate" this by
simply having a regular global variable and ensuring that all processes set the value of the
global variable on start up. Alternatively, <tt>distributed_glshared</tt> can be used as well.

Engine / core's sync_now() operation is not yet implemented.
*/

