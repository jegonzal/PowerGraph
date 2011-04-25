// this file contains the additional docs.  
/**
\defgroup group_schedulers Schedulers
\defgroup util_internal Internal Utility Classes
\defgroup util GraphLab Utility Classes and Functions
\defgroup group_rpc GraphLab RPC
\defgroup random Random Number Generators
\mainpage 
   
   \section intro_sec Introduction
   
   GraphLab is a powerful new system for designing and implementing
   parallel algorithms in machine learning.  While the current targets
   multi-core shared memory parallel systems we are in the process of
   implementing a distributed version and plan to provide support for
   alternative parallel architectures include GPUs in the near future.

  The easiest way to pick up GraphLab is to code!
  Here is a  \ref pagerank_example "pagerank example" which will provide
  you with some of the high level ideas of GraphLab
  And here is a more \ref detailed_example "detailed example" which provides
  more details as well as the the supporting APIs surrounding GraphLab.
  
  
  The key pages of interest are:
  \li The \ref graphlab::graph data structure. \n
    represents a directed graph container and is used extensively throughout GraphLab. 
  \li \ref Scopes
  \li \ref Schedulers
  \li \ref shared_data
  
  GraphLab is heavily templatized and the following two structures
  help to simplify template usage.
  \li The \ref graphlab::core data structure. \n
    This provides a convenient wrapper around most of Graphlab.
  \li \ref graphlab::types \n
    This provides typedefs for all shared memory GraphLab types.

  GraphLab is additionally supported by a serialization library, a fast parallel/thread-safe
 random number generator, as well as a flexible command line parsing system.
  \li \ref Serialization
  \li \ref random
  \li \ref graphlab::command_line_options \n
       The \ref detailed_example "detailed example" provides a good example of how this is used
       
  The GraphLab library also has a collection of \ref util "parallel utility classes and functions"
  which may be useful.

*/
