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




*/
