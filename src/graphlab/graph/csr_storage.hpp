#include <boost/range.hpp>

namespace graphlab {
  /**
   * Compressed Sparse Row Representation of a Graph Structure 
   * id must be consecutive integers (or equivalent types).
   */
  template<typename IDTYPE>
  class csr_storage {
    typedef IDTYPE id_type;
    typedef size_t ptr_type;

   public:
    size_t num_entries() {
      return ids.size();
    }

    std::vector<ptr_type> get_index() {
      return index;
    }

    std::vector<id_type> get_ids() {
      return ids;
    }

    // Private Data Storage
   private:
    std::vector<ptr_type> index;
    std::vector<id_type> ids; 
  } // end of csr_storage 
} // end of graphlab
