#include <graphlab/factors/table_factor.hpp>


std::ostream& operator<<(std::ostream& out, const graphlab::variable& var) {
  // return out << "v_" << var.id
  //            << " in {0:" << var.arity-1 << "}";
  return out << var.id;
}

