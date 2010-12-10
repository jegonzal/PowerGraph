// this file contains the additional docs.  
/**
\defgroup group_schedulers Schedulers
\defgroup util_internal Internal Utility Classes


\mainpage 
   
   \section intro_sec Introduction
   
   GraphLab is a powerful new system for designing and implementing
   parallel algorithms in machine learning.  While the current targets
   multi-core shared memory parallel systems we are in the process of
   implementing a distributed version and plan to provide support for
   alternative parallel architectures include GPUs in the near future.
 
   For a more user friendly tour of GraphLab and its features visit
   the site: <a href="http://www.graphlab.ml.cmu.edu/details.html">
   http://www.graphlab.ml.cmu.edu/details.html </a>
    Those more interested in learning about the abstraction should read
   the <a href="http://www.select.cs.cmu.edu/publications/scripts/papers.cgi?Low+al:uai10graphlab">
   GraphLab paper </a>.

\endpage


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
    type can be selected through the \ref iengine::set_default_scope()
    function in the \ref engine, or the
    ref core::set_scope_type() function in the \ref core .
    
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
      <dd>This scope type guarantees safe read/write access 
      to the current vertex as well as adjacent edges. In addition,
      you can also read (but not write) data on adjacent vertices safely.
      This is implemented by acquiring a write lock on the current
      vertex, and a read lock on all adjacent vertices. 
    \li <b> Vertex Consistency </b>
      <dd>This scope type guarantees only safe read/write access 
      to the current vertex. This is implemented by acquiring a write
      lock on the current vertex. 
\endpage


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
\li <b> \ref core::sched_options() </b>
If the \ref core object is used, the scheduler options can be set using
the core::sched_options() function by calling
<tt> core.sched_options().add_option("optionname", VALUE) </tt>
where VALUE can be any type. (sched_options will cast to the right type as
needed automatically).

\section sec_static_schedulers Static Schedulers

 The static schedulers ignore
  add_task calls, and typically repeat a predefined
  sequence until termination.

\subsection sec_colored_scheduler Colored Scheduler ["colored_scheduler"]

      The colored scheduler supports only a single update function and
      applies it in sweeps over the graph.  The colored reads the
      color of each vertex from the \ref graph::color()
      function.  The scheduler then executes all vertices of the same
      color before proceeding to the next color.  Within a color,
      processors may execute update functions on multiple vertices
      simultaneously.  However, the colored scheduler ensures that at
      no point are vertices of different colors executed
      simultaneously. 

      The color scheduler can be used to obtain a Gauss-Seidel
      parallel execution of a sequential algorithm by first computing
      a coloring of the underlying model
      using \ref graph::compute_coloring().  As long as the
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
      effect on the scheduling, it is used here by the colored_scheduler
      to provide a default update function to use.

\subsection sec_round_robin Round Robin Scheduler ["round_robin"]

      This scheduler executes repeated cycles through all vertices.
      The add_task() functions can be used to associate
      different update functions with each vertex.
    
\subsection sec_set_scheduler set_scheduler ["set"]
    This is experimental and should not be used.
    It is not part of the official release.


\section sec_dynamic_schedulers Dynamic Schedulers 
  
    The dynamic schedulers rely on the \ref icallback in the
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
         Sets the partition method to use. See \ref partition_method
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

\endpage


*/