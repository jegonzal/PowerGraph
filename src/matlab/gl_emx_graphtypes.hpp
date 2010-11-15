#ifndef GL_EMX_GRAPHTYPES_HPP
#define GL_EMX_GRAPHTYPES_HPP
#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/typeof/std/utility.hpp>
#include <boost/function.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>

#include "datatype_identifier.h"

#define F1(F) boost::function<F>::arg1_type
#define F2(F) boost::function<F>::arg2_type


typedef BOOST_TYPEOF(datatype_identifier) datatype_identifier_fn_type;

typedef F1(datatype_identifier_fn_type) gl_emx_vertextype;
typedef F2(datatype_identifier_fn_type) gl_emx_edgetype;

#undef F1
#undef F2
#endif