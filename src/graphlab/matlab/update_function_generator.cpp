#include <mex.h>
#include <boost/unordered_map.hpp>
#include <boost/preprocessor.hpp>
#include "gl_emx_graphtypes.hpp"
#include "update_function_generator.hpp"
//#include "update_function_array.hpp"
/*********************************************************************
 * The __matlab_config.hpp defines an array __UPDATE_FUNCTIONS__ 
 * which lists all the matlab update functions available.
 * However, I need to generate my own graphlab update functions 
 * which call the matlab update functions. So, for each matlab update 
 * function, I generate a graphlab update function with the same name but
 * prefixed with __gl__
 *********************************************************************/

// Some useful macros to extract elements from the __UPDATE_FUNCTIONS__ array 
#define GET_UPDATE_FUNCTION_N(N) BOOST_PP_ARRAY_ELEM(N, __UPDATE_FUNCTIONS__)
#define GET_NUM_UPDATE_FUNCTIONS BOOST_PP_ARRAY_SIZE(__UPDATE_FUNCTIONS__)
#define GET_UPDATE_FUNCTION_NAME_N(N) BOOST_PP_STRINGIZE(GET_UPDATE_FUNCTION_N(N))
// To get the GraphLab versions of the update functions
#define GET_GL_UPDATE_FUNCTION_N(N) BOOST_PP_CAT( __gl__ , BOOST_PP_ARRAY_ELEM(N, __UPDATE_FUNCTIONS__))
#define GET_GL_UPDATE_FUNCTION_NAME_N(N) BOOST_PP_STRINGIZE(GET_GL_UPDATE_FUNCTION_N(N))


#define GEN_UPDATE_FUNCTION(Z,N,_) void GET_GL_UPDATE_FUNCTION_N(N)     \
                                  (gl_types::iscope& scope,               \
                                   gl_types::icallback& scheduler,         \
                                   gl_types::ishared_data* shared_data) {  \
  exec_update_function< GET_UPDATE_FUNCTION_N(N) >(scope,scheduler,shared_data); \
}
  
BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, GEN_UPDATE_FUNCTION, _)
  
/*********************************************************************
 * Registers all the graphlab update functions created in a map
 * so I get a string->update_function mapping.
 *********************************************************************/
#define INSERT_ELEM_INTO_MAP(Z,N, _) \
  update_function_map[GET_GL_UPDATE_FUNCTION_NAME_N(N)] = GET_GL_UPDATE_FUNCTION_N(N); \
  mexPrintf("Registering update function: " GET_GL_UPDATE_FUNCTION_NAME_N(N) "\n");
update_function_map_type update_function_map;

void register_all_matlab_update_functions() {
  BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, INSERT_ELEM_INTO_MAP, _)
}

