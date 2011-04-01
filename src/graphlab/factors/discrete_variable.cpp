#include <graphlab/factors/discrete_variable.hpp>
std::ostream& operator<<(std::ostream& out, 
                         const graphlab::discrete_variable& var) {
  // return out << "v_" << var.id()
  //            << " in {0:" << var.size()-1 << "}";
  return out << var.id();
}


