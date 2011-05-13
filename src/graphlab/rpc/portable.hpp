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

