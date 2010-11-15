#include <boost/unordered_map.hpp>
#include <boost/preprocessor.hpp>
#include "graphlab_mex.hpp"
#include "__matlab_config.hpp"

/*********************************************************************
 * The __matlab_config.hpp defines an array UPDATE_FUNCTIONS 
 * which lists all the matlab update functions available.
 * However, I need to generate my own graphlab update functions 
 * which call the matlab update functions. So, for each matlab update 
 * function, I generate a graphlab update function with the same name but
 * prefixed with __gl__
 *********************************************************************/

// Some useful macros to extract elements from the UPDATE_FUNCTIONS array 
#define GET_UPDATE_FUNCTION_N(N) BOOST_PP_ARRAY_ELEM(N, UPDATE_FUNCTIONS)
#define GET_NUM_UPDATE_FUNCTIONS BOOST_PP_ARRAY_SIZE(UPDATE_FUNCTIONS)
#define GET_UPDATE_FUNCTION_NAME_N(N) BOOST_PP_STRINGIZE(GET_UPDATE_FUNCTION_N(N))
// To get the GraphLab versions of the update functions
#define GET_GL_UPDATE_FUNCTION_N(N) BOOST_CAT( __gl__ , BOOST_PP_ARRAY_ELEM(N, UPDATE_FUNCTIONS))
#define GET_GL_UPDATE_FUNCTION_NAME_N(N) BOOST_PP_STRINGIZE(GET_GL_UPDATE_FUNCTION_N(N))



/*********************************************************************
 * Generate all the update functions
 * To simplify the amount of code I need to generate using boost PP,
 * I simply use the boost preprocessor create a function which calls
 * (and instantiates) a generic templateized update function caller.
 *********************************************************************/

template <matlab_update_function mlupdate>
void exec_update_function(mrf_gl::iscope& scope, 
                          mrf_gl::icallback& scheduler,
                          mrf_gl::ishared_data* shared_data) {
  std::cout << "!" << std::endl;
}

#define GEN_UPDATE_FUNCTION(Z,N,_) void GET_GL_UPDATE_FUNCTION_N(N)     \
                                  (mrf_gl::iscope& scope,               \
                                  mrf_gl::icallback& scheduler,         \
                                  mrf_gl::ishared_data* shared_data) {  \
  exec_update_function< GET_UPDATE_FUNCTION_N(N) >(scope,scheduler,shared_data); \
}
  
BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, GEN_UPDATE_FUNCTION, _)
  
/*********************************************************************
 * Registers all the graphlab update functions created in a map
 * so I get a string->update_function mapping.
 *********************************************************************/
#define INSERT_ELEM_INTO_MAP(Z,N, _) update_function_map[GET_GL_UPDATE_FUNCTION_NAME_N(N)] = GET_GL_UPDATE_FUNCTION_N(N) ;
boost::unordered_map<std::string, update_task<matlab_graph>::update_function_type> update_function_map;
void register_all_matlab_update_functions() {
  BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, INSERT_ELEM_INTO_MAP, _)
}


