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

// this file contains the additional docs.  
/**
\defgroup group_schedulers Schedulers
\defgroup util_internal Internal Utility Classes
\defgroup util GraphLab Utility Classes and Functions
\defgroup rpc GraphLab RPC
\defgroup rpc_internal GraphLab RPC Internal
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

  Additionally, GraphLab also has a nice RPC implementation (which we use
  in the distributed implementation which is not yet available)
  \li \ref RPC "GraphLab RPC"
*/
