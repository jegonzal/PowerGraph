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

#include "get_edge_data.h"
#include "get_vertex_data.h"

#define F1(F) boost::function<F>::arg1_type
#define F2(F) boost::function<F>::arg2_type
#define F3(F) boost::function<F>::arg3_type
#define FRESULT(F) boost::function<F>::result_type
#define FARITY(F) boost::function<F>::arity

namespace elicit_graph_type_detail {

template <typename F>
struct has_2_args {
  typedef typename boost::mpl::bool_<FARITY(F) == 2 >::type type;
};


template <typename F, size_t nargs>
struct get_args{
 typedef typename FRESULT(F) result_type;
 typedef typename F3(F) arg3_type;
};

// if 0 args. then arg3 is void
template <typename F>
struct get_args<F, 0>{
 typedef typename FRESULT(F) result_type;
 typedef void arg3_type;
};

// if 1 args. then arg3 is void
template <typename F>
struct get_args<F, 1>{
 typedef typename FRESULT(F) result_type;
 typedef void arg3_type;
};

// if 2 args. then arg3 is void
template <typename F>
struct get_args<F, 2>{
 typedef typename FRESULT(F) result_type;
 typedef void arg3_type;
};



template <typename F>
struct get_data_ret_type {
  typedef typename boost::mpl::if_< typename has_2_args<F>::type,
               typename get_args<F,FARITY(F)>::result_type,
               typename boost::remove_pointer<typename get_args<F,FARITY(F)>::arg3_type>::type
               >::type type;
};

}

typedef BOOST_TYPEOF(get_vertex_data) get_vertex_data_fn_type;
typedef BOOST_TYPEOF(get_edge_data) get_edge_data_fn_type;

typedef elicit_graph_type_detail::get_data_ret_type<get_vertex_data_fn_type>::type gl_emx_vertextype;
typedef elicit_graph_type_detail::get_data_ret_type<get_edge_data_fn_type>::type gl_emx_edgetype;
#endif