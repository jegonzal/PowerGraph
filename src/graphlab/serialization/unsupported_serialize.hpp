/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
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

