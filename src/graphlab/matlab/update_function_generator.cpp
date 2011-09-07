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
                                   gl_types::icallback& scheduler) {  \
  exec_update_function< GET_UPDATE_FUNCTION_N(N) >(scope,scheduler); \
}
  
BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, GEN_UPDATE_FUNCTION, _)
  
/*********************************************************************
 * Registers all the graphlab update functions created in a map
 * so I get a string->update_function mapping.
 *********************************************************************/
#define INSERT_ELEM_INTO_MAP(Z,N, _) \
  update_function_map[GET_GL_UPDATE_FUNCTION_NAME_N(N)] = GET_GL_UPDATE_FUNCTION_N(N); \
  std::cout << "Registering update function: " << GET_GL_UPDATE_FUNCTION_NAME_N(N) << std::endl;
update_function_map_type update_function_map;

void register_all_matlab_update_functions() {
  BOOST_PP_REPEAT(GET_NUM_UPDATE_FUNCTIONS, INSERT_ELEM_INTO_MAP, _)
}

