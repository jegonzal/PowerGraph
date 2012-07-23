/*  
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





#ifndef GRAPHLAB_CHANDY_MISRA_INTERFACE_HPP
#define GRAPHLAB_CHANDY_MISRA_INTERFACE_HPP

namespace graphlab {
  template <typename GraphType>
  class chandy_misra_interface {
  public:
    typedef typename GraphType::lvid_type lvid_type;
    virtual ~chandy_misra_interface() { }
    virtual size_t num_clean_forks() const = 0;
    virtual void make_philosopher_hungry_per_replica(lvid_type p_id) = 0; 
    virtual void philosopher_stops_eating_per_replica(lvid_type p_id) = 0;
    
  };
} // end of GraphLab namespace
#endif
