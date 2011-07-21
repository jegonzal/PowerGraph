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


#ifndef GRAPHLAB_EXECUTION_STATUS_HPP
#define GRAPHLAB_EXECUTION_STATUS_HPP

namespace graphlab {

  /**
   * \brief the reasons for execution completion.
   *
   * Because there are several reasons why the graphlab engine might
   * terminate the exec_status value is returned from the start
   * function after completing execution. 
   *
   */
  struct execution_status {
    enum status_enum {
      EXEC_UNSET,  /** The default termination reason */      
      EXEC_TASK_DEPLETION, /**<Execution completed successfully due to
                              task depletion */      
      EXEC_TERM_FUNCTION,  /**< Execution completed successfully due to
                              termination function. */      
      EXEC_TIMEOUT,       /**< The execution completed after timing
                             out */
      EXEC_TASK_BUDGET_EXCEEDED, /**< The execution completed because
                                    the maximum number of tasks was
                                    exceeded */
      
      EXEC_FORCED_ABORT,     /**< the engine was stopped by calling force
                                abort */
      
      EXEC_EXCEPTION        /**< the engine was stopped by an exception */
    }; // end of enum
    
    // Convenience function.
    static std::string to_string(status_enum es) {
      switch(es) {
      case EXEC_UNSET: return "engine not run!";
      case EXEC_FORCED_ABORT: return "forced abort";
      case EXEC_TASK_BUDGET_EXCEEDED: return "budget exceed";
      case EXEC_TERM_FUNCTION: return "termination function";
      case EXEC_TASK_DEPLETION: return "task depletion (natural)";
      case EXEC_TIMEOUT: return "timeout";
      case EXEC_EXCEPTION: return "exception";
      };
      return "unknown";
    } // end of to_string
  };



}; // end of namespace graphlab
#endif



