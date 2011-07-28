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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef GRAPHLAB_IAGGREGATOR_HPP
#define	GRAPHLAB_IAGGREGATOR_HPP


#include <graphlab/logger/logger.hpp>
#include <graphlab/scope/iscope.hpp>




namespace graphlab {

  /**
   * \brief  This class is the base type of an aggregation (sync)
   * operation. 
   *
   * Each thread is assumed to have a a single iaggregator which it
   * uses to assemble the the partial sum.  This class is not thread
   * safe and therefore proper locking should be used.
   *
   * The iaggregator maintains an internal partial sum.
   *
   */
  template<typename Graph> 
  class iaggregator {
  public:
    typedef Graph graph_type;
    typedef iscope<graph_type> iscope_type;

    virtual ~iaggregator();


    /**
     * Make a virtual copy of the aggregator
     */
    virtual iaggregator* clone() = 0;

    /**
     * Clear the current aggregator
     */
    virtual void clear() = 0;

    /**
     * Add the scope to the current partial sum
     */
    virtual void operator+=(const iscope_type& scope) = 0;

    /**
     * Add another partial sum to this partial sum.
     */
    virtual void operator+=(const iaggregator& other) = 0;

    /**
     * Apply this partial sum to a global shared object.
     */
    virtual void apply() = 0;
  }; // end of iaggregator




  
  
}; // end of namespace graphlab


#endif	/* GRAPHLAB_IAGGREGATOR_HPP */

