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

#ifndef GRAPHLAB_MEMORY_INFO_HPP
#define GRAPHLAB_MEMORY_INFO_HPP

namespace graphlab {
  /**
   * \internal \brief Memory info namespace contains functions used to
   * compute memory usage.
   *
   * The memory info functions require TCMalloc to actually compute
   * memory usage values. If TCMalloc is not present then calls to
   * memory info will generate warnings and return the default value.
   */
  namespace memory_info {

    /**
     * \interanl 
     *
     * \brief Returns whether memory info reporting is
     * available on this system (if memory_info was built with TCMalloc)
     *
     * @return if memory info is available on this system.
     */
    bool available();

    /**
     * \internal
     * 
     * \brief Estimates the total current size of the memory heap in
     * bytes. If memory info is not available then 0 is returned.
     *
     * @return size of heap in bytes
     */
    size_t heap_bytes();

    /**
     * \internal
     *
     * \brief Determines the total number of allocated bytes.  If
     * memory info is not available then 0 is returned.
     *
     * @return the total bytes allocated 
     */
    size_t allocated_bytes();

    /**
     * \internal
     * 
     * \brief Print a memory usage summary prefixed by the string
     * argument.
     *
     * @param [in] label the string to print before the memory usage summary.
     */
    void print_usage(const std::string& label = "");

    /**
     * \internal
     * 
     * \brief Log a memory usage summary prefixed by the string
     * argument.
     *
     * @param [in] label the string to print before the memory usage summary.
     */
    void log_usage(const std::string& label = "");
  } // end of namespace memory info
};

#endif


