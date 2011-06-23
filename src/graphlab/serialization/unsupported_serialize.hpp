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


#ifndef GRAPHLAB_UNSUPPORTED_SERIALIZE_HPP
#define GRAPHLAB_UNSUPPORTED_SERIALIZE_HPP

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {

  struct unsupported_serialize {
    void save(oarchive& archive) const {      
      ASSERT_MSG(false, "trying to serialize an unserializable object");
    }
    void load(iarchive& archive) {
      ASSERT_MSG(false, "trying to deserialize an unserializable object");
    }
  }; // end of struct
};


/**
Disables serialization of a class so that it will fault at runtime.
Must be declared in the global namespace.
*/
#define GRAPHLAB_UNSERIALIZABLE(tname) \
  BEGIN_OUT_OF_PLACE_LOAD(arc, tname, tval) \
    ASSERT_MSG(false,"trying to deserialize an unserializable object"); \
  END_OUT_OF_PLACE_LOAD()                                           \
  \
  BEGIN_OUT_OF_PLACE_SAVE(arc, tname, tval) \
    ASSERT_MSG(false,"trying to serialize an unserializable object"); \
  END_OUT_OF_PLACE_SAVE()                                           \


#endif

