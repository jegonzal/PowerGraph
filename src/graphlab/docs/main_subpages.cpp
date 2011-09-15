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
\page Scopes Scopes
    The \b scope of a vertex is the data on the vertex, the data
    on all adjacent edges as well as data on all adjacent vertices. An
    update function executed on a vertex has the ability to read and
    write to all data within the scope of the vertex.

    The ability for an update function to read or write a piece of
    data does not mean that it can do so in a \b consistent
    fashion. For instance, without adequate protection, it is
    conceivable that update functions executed on neighboring vertices
    could modify the same piece of data at the same time.

    GraphLab therefore provides the concept of a \b consistency
    class which expresses both the level of protection provided, as
    well as the amount of parallelism available. The consistency
    type can be selected through the \ref graphlab::iengine::set_default_scope()
    function in the \ref graphlab::iengine, or the
    ref graphlab::core::set_scope_type() function in the \ref graphlab::core .
    
    There are three basic types of scopes (as well a number of
    undocumented ones which you probably should not use).
    The scopes are currently implemented through the use of a
    read-write lock on each vertex. Different types of scopes
    then make use of different locking schemes.

    \li <b> Full Consistency </b>
      This scope type is the safest and guarantees safe access 
      to the current vertex, adjacent edges, and adjacent vertices.
      This is implemented by acquiring a write lock on all vertices
      in the scope.
     \li <b> Edge Consistency </b>
      This scope type guarantees safe read/write access 
      to the current vertex as well as adjacent edges. In addition,
      you can also read (but not write) data on adjacent vertices safely.
      This is implemented by acquiring a write lock on the current
      vertex, and a read lock on all adjacent vertices. 
    \li <b> Vertex Consistency </b>
      This scope type guarantees only safe read/write access 
      to the current vertex. This is implemented by acquiring a write
      lock on the current vertex. 



\page Schedulers Schedulers
The GraphLab framework provides a collection of built-in 
basic schedules as well as some more exotic specialized dynamic schedules.  
Here, we provide a brief description of each scheduler as well as
the options available for each scheduler.

There are many ways to provide options to the scheduler:
\li <b> Command Line Options </b>
If the graphlab command_line_options parser is used, many of these options
can be passed through the command line using the option \n
<tt> --scheduler schedtype[key=value,key=value...] </tt> \n
Help about the schedulers can also be obtained from the <tt> --schedhelp </tt>
option.
\li <b> \ref graphlab::core::sched_options() </b>
If the \ref graphlab::core object is used, the scheduler options can be set using
the core::sched_options() function by calling
<tt> core.sched_options().add_option("optionname", VALUE) </tt>
where VALUE can be any type. (sched_options will cast to the right type as
needed automatically).

\section sec_static_schedulers Static Schedulers

 The static schedulers ignore
  add_task calls, and typically repeat a predefined
  sequence until termination.

\subsection sec_chromatic_scheduler Chromatic Scheduler ["chromatic_scheduler"]

      The chromatic scheduler supports only a single update function and
      applies it in sweeps over the graph.  The chromatic reads the
      color of each vertex from the \ref graphlab::graph::color()
      function.  The scheduler then executes all vertices of the same
      color before proceeding to the next color.  Within a color,
      processors may execute update functions on multiple vertices
      simultaneously.  However, the chromatic scheduler ensures that at
      no point are vertices of different colors executed
      simultaneously. 

      The color scheduler can be used to obtain a Gauss-Seidel
      parallel execution of a sequential algorithm by first computing
      a coloring of the underlying model
      using \ref graphlab::graph::compute_coloring().  As long as the
      update function does not \b modify the data on neighboring
      vertices, the parallel execution will be identical to all
      sequential executions in which the vertices are updated in order
      of color.

      Options:
      \li \b "max_iterations" [integer, default = 0]
           Sets the number of iterations. 
      \li \b "update_function" [update_function_type, default = set on add_task]
      Set the update function to call on the vertices. Since this takes a 
      function, this parameter cannot be filled on the command line, but must
      be filled in program code. If this parameter is not provided,
      the update function will be set to the update function provided on 
      the most recent call to add_task(). Even though add_task() has no
      effect on the scheduling, it is used here by the chromatic_scheduler
      to provide a default update function to use.

\subsection sec_round_robin Round Robin Scheduler ["round_robin"]

      This scheduler executes repeated cycles through all vertices.
      The add_task() functions can be used to associate
      different update functions with each vertex.
    
\subsection sec_set_scheduler set_scheduler ["set"]
    This is experimental and should not be used.
    It is not part of the official release.


\section sec_dynamic_schedulers Dynamic Schedulers 
  
    The dynamic schedulers rely on the \ref graphlab::icallback in the
    update functions to receive new tasks (and potentially task
    priorities).  These tasks are then incorporated into the execution
    schedule.  These schedulers are called dynamic schedulers because
    they rely on the computation in the update functions to determine
    what to execute and when.  

    All of the dynamic schedulers have an internal task de-duplication
    mechanism that prevents the same task from occurring more than
    once in the queue.  Schedulers that use priorities typically take
    the maximum priority of all duplicate entries and assign that to a
    single instance in the queue.  There are quite a few dynamic
    schedulers in GraphLab and below is a list of several of the most
    popular:

\subsection sec_fifo_scheduler FiFo Scheduler ["fifo"]  

The fifo scheduler executes tasks in the classical first in
      first out order.  When update functions generate tasks they are
      added to the back of the fifo queue.  Because there is a single
      central queue the fifo scheduler can become a synchronizing
      bottleneck for algorithms with relatively light update
      functions.  

\subsection sec_multiqueue_fifo_scheduler Multiqueue FiFo Scheduler ["multiqueue_fifo"] 
      The Multiqueue fifo scheduler is like the fifo scheduler but
      instead of using a single queue multiple queues (2 x ncpus) are
      used and tasks are added to queues using a randomized balancing
      strategy.  Each processor only draws from its own pair of queue.

\subsection sec_priority_scheduler Priority Scheduler ["priority"]

      The priority scheduler maintains a single priority scheduling
      queue.  The task with highest priority is always executed next.
      If add_task() is invoked on an already present task
      the existing task's priority is set to the max of the two tasks.


\subsection sec_multiqueue_priority_scheduler Multiqueue Priority Scheduler ["multiqueue_priority"]

      Same as the priority scheduler except multiple queues are
      maintained.


\subsection sec_clustered_priority_scheduler Clustered Priority Scheduler ["clustered_priority"]

      The clustered priority schedule maintains a priority queue over
      blocks of vertices.  This schedule begins by partitioning the
      graph using the partitioning method
      into blocks which each contain vert_per_part vertices.
      The priorities are maintained over blocks but within each block
      a sweep schedule is run.

      Options:
     \li \b "partition_method" [string: metis/random/bfs, default=metis]
         Sets the partition method to use. See \ref graphlab::partition_method
         for details about each partitioning method
     \li "vertices_per_partition" [integer, default = 100]
         Number of vertices in each partition

\subsection sec_sweep_schedueler Sweep Scheduler ["sweep"] 
      This scheduler loops over vertices executing a task if one is
      associated with the vertex.  Each vertex maintains its own local
      queue.  This scheduler has the least possible overhead.

     Options:
     \li \b "ordering" [string: linear/permute, default=linear]
     Sets the ordering to use. If "linear", the scheduler will
     sweep over all the vertices in order of the vertex id. If "permute",
     the scheduler will sweep over the vertices in a permuted order.

\section sec_specialized_schedulers Specialized Schedulers 
  \subsection sec_splash_scheduler Splash Scheduler ["splash"]
  
  We currently only provide the Splash scheduler which grows small
    spanning trees for each cpu and then sequential executes the
    update function in a forward backward fashion.  This scheduler is
    loosely based on the Splash BP algorithm by Gonzalez et
    al. AISTATS 2009. 
   
  Options:
    \li \b "splash_size" [integer, default = 100]
    The number of vertices in a splash.
    \li \b "update_function" = [update_function_type,default = set on add_task_to_all]
    The update function to use for the Splash update. Since this takes a 
      function, this parameter cannot be filled on the command line, but must
      be filled in program code. If this parameter is not provided,
      the update function will be set to the update function provided on 
      the most recent call to add_task_to_all(). 





\page Serialization Serialization
We have a custom serialization scheme which is designed for performance rather than 
compatibility. It does not perform type checking, It does not perform pointer tracking, 
and has only limited support across platforms. It has been tested, and should be compatible 
across x86 platforms.

There are two serialization classes graphlab::oarchive and graphlab::iarchive. The former 
does output, while the latter does input. To include all serialization headers,
#include <graphlab/serialization/serialization_includes.hpp>.

\section sec_basic_serialize Basic serialize/deserialize

To serialize data to disk, you just create an output archive, and associate it with a file stream.

\code
std::ofstream fout("file.bin", std::fstream::binary);
graphlab::oarchive oarc(fout);
\endcode

The stream operators are then used to write data into the archive.

\code
int i = 10;
double j = 20;
std::vector v(10,1.0);

oarc << i << j << v;
\endcode

To read back, you use the iarchive with an input file stream, and read 
back the variables in the same order:

\code
std::ifstream fin("file.bin", std::fstream::binary);
graphlab::iarchive iarc(fout);
  
int i;
double j;
std::vector v;

iarc >> i >> j >> v;
\endcode

All basic datatypes are supported. 
The four STL containers std::map, std::set, std::list and std::vector are supported 
as long as the contained type can be serialized. That is to say, it will correctly 
serialize a vector of a set of integers.

\section sec_serialize_user User Structs and Classes

To serialize a struct/class, all you need to do is to define a public load/save function. For instance:

\code
class TestClass{
 public:
  int i, j;
  std::vector<int> k;

  void save(oarchive& oarc) const {
    oarc << i << j << k;
  }

  void load(iarchive& iarc) {
    iarc >> i >> j >> k;
  }
};
\endcode

After which, the standard stream operators as described in the previous section 
will work fine. STL containers of TestClass will work as well.

\section sec_serialize_pod POD (Plain Old Data) types
POD types are data types which occupy a contiguous region in memory. For instance, basic types 
(double, int, etc), or structs which contains only basic types. Such types can be copied or
replicated using a simple mem-copy operation and is a great candidate for acceleration during
serialization / deserialization.

Unfortunately, there is not a simple way to detect or test if any given type is POD or not.
C++ TR1 defines an is_pod feature, but this feature is yet implemented in some compilers.
We therefore defined our own gl_is_pod<...> feature which allows the user to explicitly
express that a particular type is a POD type. gl_is_pod can be further extended in the future
when compilers provide better support for std::is_pod.

To use gl_is_pod, we consider the following Coordinate struct. 
\code
struct Coordinate{
  int x, y, z;
};
\endcode

This struct can be defined to be a POD type using an accelerated serializer by defining:

\code
namespace graphlab {
  struct gl_is_pod<Coordinate> {
    BOOST_STATIC_CONSTANT(bool, value = true);
  };
}
\endcode

Now, Coordinate variables, or even vector<Coordinate> variables will serialize/deserialize faster,
making use of direct mem-copies. 

A caveat of this fast serialization mechanism is that the resulting archive may not be cross platform
since it forces a particular length for the integers inside the struct. There may also be issues 
using the same archive between compilers or different alignment options. 


\section sec_serialization_any Any
graphlab::any is a variant type derived from Boost Any, but extended to support the GraphLab
serialization system. See graphlab::any for details.






\page shared_data Shared Data

The shared data system provides controlled thread-safe access to  
global variables. They are mainly used to provide two capabilities:

\li <b> Globally Shared Variables </b>
    Allows all GraphLab computation to have controlled access to "global variables". 
    Indeed in the shared memory parallelism case, this can be said to be redundant. 
    However, the abstraction provided here extends to the distributed memory case.
\li <b> Sync </b>
  A fold/reduction framework which performs background accumulation of all vertex data.

The \ref detailed_example "detailed example" provides a good description of how this is used. 
  
\section shared_Data_gsv Globally Shared Variables

All shared data must be global variables which are declared using the syntax. 

\code
gl::glshared<T> var;
\endcode

, where T is the data type of variable. Arbitrary data types 
can be used and is guaranteed to be thread-safe. The current implementation operates 
by keeping two reference counted copies of the data and so will consume twice the 
memory capacity of using a regular variable.

To read the value of the variable by value:

\code
T val = var.get_val();
\endcode


To read the value of the variable by reference:

\code
boost::shared_ptr<T> ptr = var.get_ptr();
\endcode 

The get_ptr() function returns a pointer to the data. 
The pointer is in the form of a boost::shared_ptr which can be dereferenced 
like a usual pointer using the dereference operator (*). While the pointer 
is still in scope, the values read from the pointer is guaranteed to not change. 
Under particular conditions, holding on to the shared pointer could prevent writes 
to the variable from progressing. The pointer should therefore be released as soon 
as possible by either letting it go out of scope, or by calling ptr.release() explicitly.

The variable is modified using the set_val() function, which is quite self-explanatory.

\code
var.set_val(T newval);
\endcode

The shared data type also provides Atomic operations.

\subsection sec_gsv_atomics Atomics

In addition to the regular get/set functions, the shared variable also provides two basic atomic operations.

\code
T oldval = glshared<T>::exchange(T newval)
\endcode

Writes the "newvalue" into the variable while returning a copy of the previous value.

\code
void glshared<T>::apply(apply_function_type fun, const any& closure)
\endcode

Calls the function fun with a reference to the value of the shared variable. The type of the function is

\code
void(*apply_function_type)(any& current_value,
                           const any& closure);
\endcode
It is unfortunate, but due to current design limitations, the current_value must be passed 
to the apply function as an any type. However, the type of the value within the any container 
is the type of original glshared variable. For instance, if shared variable is declared using:

\code
      gl::glshared<size_t> shared_counts;
\endcode

And I write an apply function which adds the closure value to the shared_counts, 
the resultant code will be:

\code
      void add_counts(any& current_value, const any& increment) {
        current_value.as<size_t>() += increment.as<size_t>();
      }
\endcode

Observe that the current_value can be accessed (and written to) through the 
\ref graphlab::any::as() function and the type to be passed to the as() function 
matches the type of the shared variable. The apply function is allowed to make 
changes to the current_value, and future reads from this variable will return the new value.

\section sec_shared_data_sync Sync

Sync performs a function accumulation (also called fold/reduce) across all the 
vertices in the graph and writes the result to an associated entry in the shared 
data table. A Sync operation is defined by a pair of functions. A Sync function 
and a Apply function. This concept is best explained with an example. For instance, 
consider a graph with an double on each vertex.

\code
typedef graphlab::graph<double, double> graph_type;
typedef graphlab::types<graph_type> gl;
\endcode

and I would like to compute the L2 Norm of all the integers i.e. 
(square root of the sum of squares). I would define a shared variable used to
hold the final result.

\code
      gl::glshared<double> l2norm_value;
\endcode

As well as the following sync and apply functions:

\code
void squared_sum_sync(gl::iscope& scope, 
                 graphlab::any& accumulator) { 

  accumulator.as<double>() += 
         scope.const_vertex_data() * scope.const_vertex_data();
}

void square_root_apply(graphlab::any& current_data,  
                 const graphlab::any& new_data) {

  current_data.as<double>() = sqrt(new_data.as<double>());
}
\endcode

The variable is associated with the sync/apply functions using:

\code
core.set_sync(l2norm_value,       // shared variable
              squared_sum_sync,   // sync function
              square_root_apply,  // apply function
              double(0.0),        // initial sync value
              100);               // sync frequency

\endcode

When evaluated, the sync function (squared_sum_sync) will be called on the 
scope of each vertex in turn. The accumulator is passed from one function 
call to the next function call, accumulating the values as it goes along. 
The squared_sum function is passed the index of the associated entry in the 
shared data table, as well as a reference to the table and the scope of 
current vertex being evaluated.

When all the vertices are done, the apply function (square_root_apply) is 
evaluated on the result of the accumulation. The new_data variable 
contains the final value of the accumulator during the accumulation stage. 
At this point current_data is a reference to the shared variable l2norm_value. 
That is to say: the value of current_data is equivalent to the value of 
<code>l2norm_value.get_val()</code>. 
Modifications to current_data will be reflected in future accesses 
to l2norm_value.



A sync is created using the \ref graphlab::iengine::set_sync() member function of the engine
or the core \ref graphlab::core::set_sync()
(the core simply forwards the call to the engine it contains).




\subsection sec_common_sync_applys Common Syncs and Applyies

A collection of commonly used Syncs and Applies are provided in 
gl::glshared_sync_ops and gl::glshared_apply_ops.

To use the glshared_sync_ops, you must provide an accessor function of the form 
AccumulationType accessor(const vertex_data& v);

Then you can make use of the following sync functions:

\li \ref graphlab::glshared_sync_ops::sum
Adds the value on each vertex using the accessor function to read the vertex data.
\li \ref graphlab::glshared_sync_ops::l1_sum
Adds the absolute value of each vertex using the accessor function to read the vertex data.
\li \ref graphlab::glshared_sync_ops::l2_sum
Adds the squared value of each vertex using the accessor function to read the vertex data.
\li \ref graphlab::glshared_sync_ops::max
Computes the maximum value of all the vertices using the accessor function to read the vertex data.
A collection of apply function are also provided.

\li \ref graphlab::glshared_apply_ops::identity
Stores the accumulated value to the shared variable.
\li \ref graphlab::glshared_apply_ops::identity_print
Stores the accumulated value to the shared data entry and writes it to screen
\li \ref graphlab::glshared_apply_ops::increment
Adds the final accumulated value to the associated shared variable.
\li \ref graphlab::glshared_apply_ops::decrement
Subtracts the final accumulated value from the associated shared variable.
\li \ref graphlab::glshared_apply_ops::sqrt
Stores the square root of the accumulated value to the shared variable.
For instance, the following code will create a sync which behaves exactly the same way as the example above:

\code
double accessor(const double& v) {
  return v;
}

shared_data.set_sync(l2norm_value, 
                     gl::glshared_sync_ops::l2sum<double, accessor>,
                     gl::glshared_apply_ops::sqrt<double>,
                     double(0.0), 
                     100);
\endcode

\subsection sec_sync_merge Merge 

The Sync operation as described cannot be parallelized.
To permit parallelism, an additional \b merge function which combines partial synced results
must be defined. 

If a merge function is defined, the set of vertices will be partitioned into a collection
of disjoint sets. The sync function is then performed on each set in parallel, producing
a number of partial results. These partial results are then combined with the merge function.
Finally, the apply function is executed on the final result and written to the shared variable.

For instance, in the earlier sync example where the L2 Norm of the graph vertex data is computed,
intermediate results can be combined using:

\code
static void sum_merge(any& dest,
                const any& src) {
  dest.as<double>() += src.as<double>();
}
\endcode

The set Sync call must similarly be updated

\code
core.set_sync(l2norm_value,       // shared variable
              squared_sum_sync,   // sync function
              square_root_apply,  // apply function
              double(0.0),        // initial sync value
              100,                // sync frequency
              sum_merge);         // merge function
\endcode

Similarly, a collection of common Merges are provided in gl::glshared_merge_ops

\li \ref graphlab::glshared_merge_ops::sum
Sums the intermediate results.
\li \ref graphlab::glshared_merge_ops::l1_sum
Sums the absolute value of the intermediate results.
\li \ref graphlab::glshared_merge_ops::l2_sum
Sums the squared value of the intermediate results.
\li \ref graphlab::glshared_merge_ops::max
Returns the max of the intermediate results

Example:
\code
shared_data.set_sync(l2norm_value, 
                     gl::glshared_sync_ops::l2sum<double, accessor>,
                     gl::glshared_apply_ops::sqrt<double>,
                     double(0.0), 
                     100,
                     gl::glshared_merge_ops::sum<double>);

\endcode

*/

/**
\page parallel_object_intricacies Parallel Object Intricacies
All the pthread wrapper objects are built to be easy to use for the "typical" case.
i.e. 
\code
  mutex m;
  m.lock()
  ...
  m.unlock()
\endcode
Will work as expected.

What is a little more unusual is that
\code
  std::vector<mutex> m;
  m.resize(100);
\endcode
Will also correctly produce a collection of 100 distinct mutexes. 
In particular, the parallel objects can be stored in any container type and accessed directly.

To support this, the parallel objects implement some unusual behavior for the copy constructor
and operator=.

In particular, the copy constructor for the parallel objects (say the mutex), does not actually 
perform a copy, but simply performs regular construction. The reason for this design is that the 
objects are themselves not intended to be copyable (otherwise complicated reference counting 
is necessary); but yet again a copy constructor has to be specified for correct behavior of the mutex
in STL containers. In particular, on some STL implementations, the vector resize function
is equivalent to
\code
  std::vector<mutex> m;
  m.resize(100, mutex());
\endcode
and the copy constructor is used to construct vector entries, using the default constructed
mutex() as a template. Of course, to actually perform a copy is undesirable, yet again the copy 
constructor must be  implemented. Therefore we provide a copy constructor which simply performs
regular construction.

Next, operator= is implemented, but is an empty function. operator= is necessary due to default 
operator= instantiation in structs and classes. In particular:
\code
  struct s {
    int i;
    mutex m;
  }
\endcode
implicitly generates a function
\code
s& s::operator=(const s& other) {
  i = other.i;
  m = other.m;
}
\endcode
Of course given that the struct contains a mutex, s::operator= should not be used anyway. However,
the function may still be instantiated and we would like to avoid a compile-time error. We therefore
define a mutex::operator= which does nothing. This technically permits the following behavior:
\code
  struct protected_value {
    int i;
    mutex lock;
  };
  protected_value a, b;
  a.i = 10;
  b = a;
\endcode
In this case b.i will be assigned the value of a.i, but the mutex will not be copied, but is retained. 
Even though, this is the desired behavior here, we do not recommend making use of it. operator=
and the copy constructor for the struct protected_value should be explicitly written.

This affects the
 - \ref graphlab::rwlock
 - \ref graphlab::mutex
 - \ref graphlab::conditional
 - \ref graphlab::spinlock
 - \ref graphlab::simple_spinlock
 - \ref graphlab::semaphore
 - \ref graphlab::barrier
 - \ref graphlab::cancellable_barrier
*/

