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
\page RPC GraphLab RPC

GraphLab RPC primary design goal was to provide a convenient and easy to use
asynchronous communication system between \b identical binaries running
on different machines over a distributed network. It therefore provides 
MPI-like capabilities together with RPC functionality. The GraphLab distributed 
implementation is built on top of this RPC library.

GraphLab RPC uses extensive template meta-programming techniques to provide
an \b IDL-free (http://en.wikipedia.org/wiki/Interface_description_language) 
RPC system, allowing arbitrary functions to be called on program running on 
remote machines (Note that all machines must the running the same binary).

For instance, this is a particularly interesting example:
\code
#include <iostream>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;

int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  
  if (dc.procid() == 0 && dc.numprocs() >= 2) {
    dc.remote_call(1, printf, "%d + %f = %s\n", 1, 2.0, "three");
  }
  dc.barrier();
}
\endcode

Once the distributed_control object is created, \ref distributed_control::procid "dc.procid()"
provides the current machine number, while \ref distributed_control::numprocs "dc.numprocs()"
provide the total number of machines.

The if-condition is therefore entered by only the first machine, which
performs a remote call to the second machine (the first argument of remote_call
is the target machine ID). The second machine will then
execute the equivalent of
\code
  printf("%d + %f = %s\n", 1, 2.0, "three");
\endcode

We will discuss the different aspects of the RPC library seperately:
\li \ref Spawning \n
         Initialization and Starting a distributed program using GraphLab RPC
\li \ref Basic_RPC \n
         Basic usage of the RPC library. Calling of simple functions.
\li \ref OOP_RPC \n
         Advanced usage of the RPC library. Creating and managing 
         distributed object contexts.







\page Spawning Spawning and Initialization
Spawning is the process of starting an instance of GraphLab RPC on seperate 
machines. GraphLab RPC supports two spawning methods: MPI or rpcexec.py 
(a script in the tools/ directory). The MPI method is strongly recommended,
though it does require all machines to have shared access to a common file
system.




\section sec_spawning_mpi Spawning with MPI
GraphLab was tested with MPICH2, but should also with OpenMPI.
Refer to the documentation for MPICH2 or OpenMPI to set up MPI and make sure
that you can run the basic test MPI programs (MPICH2 comes with an mpdringtest). 

No additional configuration is necessary to spawn a GraphLab RPC program with MPI.

The GraphLab RPC program should begin with:

\code
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;

int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);

  dc_init_param param;
  // set additional param options here
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  ...
}
\endcode

In this case, init_param_from_mpi uses the MPI ring to exchange port numbers
and set up the RPC communication layer. 
See \ref graphlab::dc_init_param "dc_init_param" for details about additional
configuration options.

\section sec_spawning_rpcexec Spawning with rpcexec.py
rpcexec.py provides a simple way to run a process on a collection of machines,
using ssh to communicate between them. <tt>rpcexec.py --help</tt> provides
some basic help.

You will first need to create a host file which is simply a list of host names
and IP addresses:
\verbatim
localhost
192.168.1.5
node2
node3
localhost
192.168.1.5
node2
node3
\endverbatim

Running <tt>rpcexec.py -n [num to start] -f [hostsfile] `command`</tt> will read the first
execute the command on the first N hosts in the hostfile. For instance in this case, running
\verbatim
rpcexec.py -n 5 -f hostsfile ls
\endverbatim
will run the <tt>ls</tt> bash command twice on the localhost, and once on 
the three nodes : 192.168.1.5, node2, node3.

rpcexec.py also supports a 'screen' (GNU Screen) mode. Running
\verbatim
rpcexec.py -s lsscreen -n 3 -f hostsfile ls
\endverbatim
will create a `screen` session with 3 windows where one window ran `ls` on the
localhost, while two other windows sshed into 192.168.1.5 and <tt>node2</tt>, 
running the `ls` on each of them. 

rpcexec.py will terminate immediately after creating the screen session.
\verbatim
screen -r lsscreen
\endverbatim
will display and resume the screen session.

If rpcexec.py will be used to spawn the program, The GraphLab RPC program should 
begin with:
\code
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
using namespace graphlab;

int main(int argc, char ** argv) {
  dc_init_param param;
  // set additional param options here
  if (init_param_from_env(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  ...
}
\endcode

Since unlike MPI spawning, there is no existing channel for communicating
port information between the machines. rpcexec.py therefore uses environment
variables to pass information to the GraphLab RPC process. The following 
two environment variables are used:
\li \b SPAWNNODES A comma seperated list of hostnames participating in the distributed program
\li \b SPAWNID: The index of the current machine into the SPAWNNODES list. First machine
has an index value of 0.

A machine will listen on the port 10000 + SPAWNID.

See \ref graphlab::dc_init_param "dc_init_param" for details about additional
configuration options.

This spawning system is less flexibile due to the fixed port numbering. For instance,
a crashed process will keep the port in TIMED_WAIT for a few minutes, preventing
the next RPC process from running. This also prevents multiple different GraphLab RPC programs
from running on the same set of the machines.

The MPI spawner is therefore the recommended method for starting the RPC system.


\page Basic_RPC Basic RPC Usage
Once the distributed_control is set up, it can be used to call functions on remote machines.
For instance in the earlier example:
\code
if (dc.procid() == 0) {
  dc.remote_call(1, printf, "%d + %f = %s\n", 1, 2.0, "three");
}
\endcode
calls printf from machine 0 to machine 1 asynchronously.

In the GraphLab RPC terminology, a \b call is a one-way remote function call, while a 
\b request is a function call which has a return value. \b calls are executed 
asynchronously and returns immediately, while \b requests will wait for completion
of the function on the remote machine.

For instance in the code below, machine 1 could print 
either "hello world", or "world hello".
\code
if (dc.procid() == 0) {
  dc.remote_call(1, printf, "hello ");
  dc.remote_call(1, printf, "world ");
}
\endcode

Remote calls complete \b immediately, irregardless of how long the function
took on the other side. For instance, processor 0 will take almost no time
running through this code.
\code
if (dc.procid() == 0) {
  dc.remote_call(1, sleep, 1);
}
\endcode

However, since requests will wait for completion and send back the reply,
this could take about a second to run.
\code
if (dc.procid() == 0) {
  dc.remote_request(1, sleep, 1);
}
\endcode

All arguments and return values will be passed by value. Any argument type
or return type can be used as long as it is \ref Serialization "serializable".
    
\section sec_rpc_dc Distributed Control RPC
\copydoc graphlab::distributed_control

\section sec_rpc_collective Collective Operations
In addition to regular RPC operations, A collection of MPI-like collective
operations are also provided. A collective operation is a function which requires
all machines to call the same function before execution can proceed.

\subsection sec_rpc_collective_barrier Barrier
One of the most useful operations is graphlab::distributed_control::barrier()
The barrier() is functionally equivalent to MPI_Barrier(). It requires all machines
to hit the barrier, before execution is allowed to resume. 
For instance in the code below, while processor 0 is busy working at compute Pi, 
all other machines will pause at the barrier and wait for the processor 0 to complete
computation and hit the barrier, before execution can proceed.
\code
if (dc.procid() == 0) {
  compute Pi to 1 million digits
}
dc.barrier();
\endcode

\subsection sec_rpc_collective_fullbarrier Full Barrier
A \ref graphlab::distributed_control::full_barrier "Full Barrier" is also provided
through distributed_control::full_barrier(). A Full Barrier is like a barrier but 
guarantees that all RPC operations sent before the barrier must complete execution.

For instance in the example below,
The full barrier guarantees that the call to set_a_to_1() must complete on 
all remote machines before execution is allowed to proceed. 
All machines will therefore print '1'.
\code
int a = 0;
void set_a_to_1() { a = 1; }

... /* in main after initialization */ ...
dc.remote_call( /* another machine */, set_a_to_1);
dc.full_barrier();
std::cout << a;
\endcode

The full_barrier is about 2-3x more costly than the regular barrier and should be
used sparingly.

\subsection sec_rpc_collective_other_collectives Other Collectives
In addition to the barrier and the full barrier, operations such as broadcast,
gather, all_gather are also provided.
Note that the implementation of these operations are not particularly efficient
as compared to native MPI implementations due to simplistic algorithm choices.

\section sec_rpc_sequentialization Sequentialization
A slightly more unusual feature of the GraphLab RPC system is the ability to
enforce sequentialization of a sequence of RPC calls. This is particularly
useful for asynchronous usages of this RPC library and can simplify code in many
cases.

For instance, in the code below:
\code
int a = 0;
void set_a_to_1() { a = 1; }
void print_a() { std::cout << a; }

... /* in main after initialization */ ...
targetmachine = (dc.procid() + 1) % dc.numprocs();
dc.remote_call(targetmachine, set_a_to_1);
dc.remote_call(targetmachine, print_a);
...
\endcode
Note that due to the asynchronous nature of the remote_call, it is possible for
<tt>print_a()</tt> to complete on the target machine, before the variable <tt>a</tt> is set to 1.
Therefore, it is possible for the output to be '0'.

A possible solution as suggested before is to change the remote_calls to remote_requests.
However, requests incur a large performance penalty due to the need to wait for replies.

Alternatively, we can use the sequentialization key system:
\code
char oldkey = distributed_control::new_sequentialization_key();

dc.remote_call(targetmachine, set_a_to_1);
dc.remote_call(targetmachine, print_a);

distributed_control::set_sequentialization_key(oldkey);
\endcode
This enforces that calls/requests made while a key is set will always be processed
by the same thread in the thread pool on the target machine, ensuring sequentialization
of the <tt>set_a_to_1</tt> and the <tt>print_a</tt> call. 

The sequentialization key is unique to each \b thread (thread-local) 
so sequentialization of RPC calls in one thread will not affect RPC calls made by other
threads.

Keys need not be created by distributed_control::new_sequentialization_key().
distributed_control::set_sequentialization_key() can also be used, giving
an arbitrary non-zero value for the key. 
\code
char oldkey = distributed_control::set_sequentialization_key(123);

dc.remote_call(targetmachine, set_a_to_1);
dc.remote_call(targetmachine, print_a);

distributed_control::set_sequentialization_key(oldkey);
\endcode

Essentially all RPC calls made using the same key value will sequentialize.
However, since new_sequentialization_key() uses a very naive key selection system,
we recommend the use of set_sequentialization_key() especially in the case of
multi-threaded code.


\page OOP_RPC Distributed Objects

GraphLab provides a "distributed object" system which simplifies the process
of designing data structures which provide distributed computation and storage.

A GraphLab distributed object is an object which is instantiated at the same
time across all machines. The object internally contains a <tt>dc_dist_object</tt>
which provides RPC communication between distributed instances.

For instance, say we run the following code using two machines:
\code
... /* in main after initialization */ ...
graphlab::dht<std::string, std::string> str_map;
dc.barrier();

if (dc.procid() == 0) {
  str_map.set("hello", "world");
}
else if (dc.procid() == 1) {
  str_map.set("something", "other");
}
dc.barrier();
std::cout << str_map.get("hello").second;
std::cout << str_map.get("something").second;
\endcode
The DHT is a distributed object which provides a distributed key/value
store (a distributed "Hash Table"). Every entry is stored at a machine 
corresponding to a hash of the key value. Note that it is created at the same time
on all the machines.  The barrier() after creation ensures that the object is
instantiated properly on all machines before utilization.

Now, after initialization, the <tt>set</tt> function of the dht will internally
hash the key value and forward it to the right machine for processing. <tt>get</tt>
is similar. However, since the distributed object system operates on \b instances,
it is possible to create multiple distributed objects easily. For instance,
the following code will create 50 different distributed key/value maps. 
str_map[15] corresponds to the same DHT when accessed on any machine.
\code
graphlab::dht<std::string, std::string> str_map[50];
dc.barrier();
\endcode

\section sec_oop_rpc_usage Usage
We will demonstrate usage of the distributed object system using a simple 
distributed Hash Table example. Note that this is a \b very \b simple implementation,
and is not entirely correct since we are going to ignore thread-safety. But it is
sufficient to demonstrate the key concepts.

\code
class string_dht { 
 private:
  std::map<int, std::string> local_storage;
  mutable dc_dist_object<string_dht> rmi;
\endcode

First, each machine needs a local data storage. In this case we will simply 
use a std::map. The key object that provides distributed access is the 
<tt>dc_dist_object\<string_dht\> rmi;</tt>. This object creates a "context" 
for remote function calls, allowing the correct remote instance to be identified.

We will now look at the string_dht constructor. The rmi object constructor
requires a reference to the underlying distributed_control object, as well
as a pointer to the current instance:
\code
 public:
  string_dht(distributed_control &dc): rmi(dc, this) {  }
\endcode

Now, to demonstrate how the RMI object is used, lets see the set() function
\code
void set(int key, const std::string &newval) {  
  procid_t owningmachine = key % rmi.numprocs();
  if (owningmachine == rmi.procid()) {
    local_storage[key] = newval;
  }
\endcode
We use a simple hash function to identify where the key-value pair should be 
stored. Observe that the RMI object provides pretty much the same functionality
as the distributed_control object, having both dc_dist_object::numprocs()
and dc_dist_object::procid(). If the data is to be stored in the current machine, 
we simply store it. Otherwise we will need to send it to a remote machine for
storage. This is the interesting case:

\code
  else {
    rmi.remote_call(owningmachine,
                    &string_dht::set,
                    key,
                    newval);
  }
}
\endcode

The RMI object supports the same family of call/request operations as 
\ref sec_rpc_dc "distributed_control"
However, it will only work with <b>member function pointers</b>. For instance in this case,
we will be calling the set() member function on the matching instance of the string_dht
object on a remote machine. (Note that the & is important and necessary)

The get() function is similar. However, we will have to use remote requests.

\code
std::string get(int key) {  
  procid_t owningmachine = key % rmi.numprocs();
  if (owningmachine == rmi.procid()) {
    return local_storage[key];
  }
  else {
    return rmi.remote_request(owningmachine,
                              &string_dht::get,
                              key);
  }
}
\endcode

As stated earlier, this code should not be used as it is due to several limitations
such as the local_storage object is not thread-safe. Since incoming RPC calls are
generally multithreaded, locks are necessary. See dht.hpp for an equivalent 
"safe" example of a simple DHT.

\section sec_oop_rpc_context Context
Essentially, the dc_dist_object object supports the identical set of operations as the 
distributed_control object, but restricted to the \b context of a single object instance.

It includes all the regular call operations:
\li \b remote_call
\li \b fast_remote_call
\li \b control_call
\li \b remote_request
\li \b fast_remote_request
\li \b control_request

Additionally, this \b context is entirely independent of the distributed_control
object, permitting its own set of collective operations such as dc_dist_object::broadcast,
graphlab::dc_dist_object::barrier, graphlab::dc_dist_object::full_barrier, etc

Since these collective operations also operate entirely within the context of the object
instance, this permits the use of parallel collectives. For instance, I could have
two objects, and each object internally spawns threads to perform distributed computation;
using the RMI object to perform collective operations which are local to the object.

In particular, the graphlab::dc_dist_object::full_barrier() is worth taking note of. 
The graphlab::distributed_control::full_barrier() ensures completion of ALL RPC calls
including calls meant for distributed objects. Its barrier is therefore
\b global to the state of the program as a while.
The graphlab::dc_dist_object::full_barrier() however, only ensures completion of all RPC
calls within the object instance. Its barrier is therefore \b local to the state
of the distributed object. This allows each distributed object to run its own
full barriers without affecting other distributed objects.

\section sec_oop_rpc_notes Final Notes
Finally, note that the RMI object can ONLY call member function pointers.
It cannot call other global functions (such as printf).
The global context can be accessed through dc_dist_object::dc() which
returns the underlying distributed_control object, which can then be used
to call global functions. For instance:
\code
rmi.dc().remote_call(1, printf, "hello ");
\endcode
*/
