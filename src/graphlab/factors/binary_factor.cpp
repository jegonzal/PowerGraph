#include <graphlab/factors/binary_factor.hpp>

std::ostream& operator<<(std::ostream& out, 
                         const graphlab::binary_factor& fact) {
  out << "Binary Factor(v_" << fact.var1() << " in {1..."
      << fact.arity1() << "}, " 
      << ", v_ " << fact.var2() << " in {1..." 
      << fact.arity2() << "})" << std::endl;
  for(size_t i = 0; i < fact.arity1(); ++i) {
    for(size_t j = 0; j < fact.arity2(); ++j) {
      out << fact.logP(i,j) << " ";
    }
    out << std::endl;
  }
  return out;
} // end of operator<<
