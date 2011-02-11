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

