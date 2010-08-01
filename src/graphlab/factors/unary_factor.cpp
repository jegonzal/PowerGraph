#include <graphlab/factors/unary_factor.hpp>

std::ostream& operator<<(std::ostream& out, 
                         const graphlab::unary_factor& fact) {
  out << "Unary Factor(" << fact.arity()  << ")"
      << std::endl;
  for(size_t i = 0; i < fact.arity(); ++i) {
    out << fact.logP(i) << " ";
  }
  out << std::endl;
  return out;
} // end of operator<<
