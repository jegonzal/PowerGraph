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

#include <graphlab/util/generics/any.hpp>

namespace graphlab {


    __any_registration_map_type& __get_registration_map() {
    static __any_registration_map_type __any_registration_map;
    return __any_registration_map;
  }

  __any_placeholder* __any_placeholder::base_load(iarchive_soft_fail &arc) {
    uint64_t idload;
    arc >> idload;
    __any_registration_map_type::iterator i = __get_registration_map().find(idload);
    assert(i != __get_registration_map().end());
    return __get_registration_map()[idload](arc);
  }

  void __any_placeholder::base_save(oarchive_soft_fail &arc) const {
    arc << get_deserializer_id();
    save(arc);
  }
  

} // end of namespace graphlab


// std::ostream& operator<<(std::ostream& out, const graphlab::any& any) {
//   return any.print(out);
// } // end of operator << for any

