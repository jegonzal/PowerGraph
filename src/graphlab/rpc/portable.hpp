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


#ifndef PORTABLE_HPP
#define PORTABLE_HPP

#include <string>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
  
/**
\ingroup rpc
Defines a simple macro called PORTABLE which is used to wrap
a function name when issuing a portable call. 
*/  
#define PORTABLE(f) graphlab::portable_call<typeof(f)*>(BOOST_PP_STRINGIZE(f)) 

template <typename F>
struct portable_call{
  typedef F f_type;
  portable_call(){}
  portable_call(std::string c): fname(c) {}
  std::string fname;
};

} // makespace graphlab

/**
\ingroup rpc_internal
Make the portable call have the right return type
*/
namespace boost {
template <typename F>
struct function_traits<graphlab::portable_call<F> > {
  typedef __GLRPC_FRESULT result_type;
};


template <typename F>
struct function<graphlab::portable_call<F> > {
  typedef __GLRPC_FRESULT result_type;
};

}

#include <graphlab/rpc/function_arg_types_undef.hpp>
#endif

