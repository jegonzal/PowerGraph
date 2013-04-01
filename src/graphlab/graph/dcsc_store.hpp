#ifndef GRAPHLAB_DCSC_STORE_HPP
#define GRAPHLAB_DCSC_STORE_HPP
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <graphlab/graph/dcsc_store_iterator.hpp>
namespace graphlab {


/**
 * Implements a DCSC matrix representation.
 */
template <typename IndexType = uint32_t, typename ValueType = uint32_t>
struct dcsc_store {
  typedef IndexType index_type;
  typedef ValueType value_type;
  typedef dcsc_store_impl::entry_list<IndexType, ValueType> entry_list;
  typedef dcsc_store_impl::entry_type<IndexType, ValueType> entry_type;
  typedef dcsc_store_impl::entry_iterator<IndexType, ValueType> entry_iterator;

  std::vector<index_type> jc;       // column indices
  std::vector<size_t> cp;    // Parallel array to jc, for the column, points to
                             // the entry in IR for which the row indices begins.
  std::vector<index_type> ir;       // A list of row indices
  std::vector<value_type> value; // A list of values on each row index


  dcsc_store() {
    clear();
  }

  /**
   * Inserts an entry into the matrix .
   * Assumes that the row/col do not already exist.
   * Returns the offset in which the data was inserted.
   * The data can be retrieved using data_at_offset(..) if no other
   * insertions have been made
   */
  size_t insert(const index_type& row,
              const index_type& col,
              const value_type& val) {
    typename std::vector<index_type>::iterator jciter =
        std::lower_bound(jc.begin(), jc.end(), col);
    size_t jcindex = std::distance(jc.begin(), jciter);
    if (jciter == jc.end() || (*jciter) != col) {
      // we need to insert a column. Also update the cp
      // array
      jc.insert(jciter, col);
      cp.insert(cp.begin() + jcindex, cp[jcindex]);
    }
    size_t rows_begin = cp[jcindex];
    size_t rows_end = cp[jcindex + 1]; // one past the end
    // find the insert location. For very sparse, a linear search is faster
    size_t insertloc = rows_begin;
    while(insertloc < rows_end && ir[insertloc] < row) ++insertloc;

    ir.insert(ir.begin() + insertloc, row);
    value.insert(value.begin() + insertloc, val);

    // advance the cp indices after the insert location by 1
    for (size_t i = jcindex + 1; i < cp.size(); ++i) {
      ++cp[i];
    }
    return insertloc;
  }


  /**
   * Prints the matrix
   */
  void print(std::ostream& out) const {
    for (size_t c = 0;c < jc.size(); ++c) {
      for (size_t p = cp[c];p < cp[c + 1]; ++p) {
        out << "(" << ir[p] << ", " << jc[c] << ") = " << value[p] << "\n";
      }
    }
  }


  /**
   * Returns the offset in the value array for which this element can
   * be found. Returns size_t(-1) if not found.
   */
  size_t index(index_type row, index_type col) const {
    typename std::vector<index_type>::const_iterator jciter =
        std::lower_bound(jc.begin(), jc.end(), col);
    if (jciter == jc.end() || (*jciter) != col) return (size_t)(-1);

    size_t c = std::distance(jc.begin(), jciter);
    for (size_t p = cp[c];p < cp[c + 1]; ++p) {
      if (ir[p] == row) {
        return p;
      }
    }
    return (size_t)(-1);
  }

  /**
   * Returns true if this matrix contains an element at this given row/col.
   * Returns false otherwise.
   */
  bool contains(index_type row, index_type col) const {
    return index(row, col) != (size_t)(-1);
  }

  /**
   * Returns the value of the matrix at this given row/col.
   * Returns 0 if not found. Note that this function cannot distinguish
   * between if the actual value is 0, or if the entry is not found.
   * Use index() or contains() for that.
   */
  value_type find(index_type row, index_type col) const {
    size_t idx = index(row, col);
    if (idx == (size_t)(-1)) return value_type(0);
    else return value[idx];
  }

  /**
   * Empties the matrix
   */
  void clear() {
    jc.clear();
    cp.clear();
    ir.clear();
    value.clear();
    cp.push_back(0);
  }


  template <typename IteratorType>
  struct indirect_sort_functor {
    IteratorType arrbegin;
    indirect_sort_functor(IteratorType arrbegin) : arrbegin(arrbegin) { }
    bool operator()(size_t a, size_t b) const {
      return (*(arrbegin + a)) < (*(arrbegin + b));
    }
  };

  /**
   * Constructs a matrix from the arguments.
   */
  template <typename IndexTypeIterator, typename ValueTypeIterator>
  void construct(IndexTypeIterator rowbegin, IndexTypeIterator rowend,
                 IndexTypeIterator colbegin, IndexTypeIterator colend,
                 ValueTypeIterator valbegin, ValueTypeIterator valend) {
    assert(std::distance(rowbegin, rowend) == std::distance(colbegin, colend));
    assert(std::distance(rowbegin, rowend) == std::distance(valbegin, valend));
    size_t numvals = std::distance(rowbegin, rowend);

    clear();
    if (numvals == 0) return;

    // generate a sorting permutation of the columns
    std::vector<index_type> indices(std::distance(colbegin, colend));
    for (size_t i = 0;i < indices.size(); ++i) indices[i] = i;
    std::sort(indices.begin(), indices.end(),
              indirect_sort_functor<IndexTypeIterator>(colbegin));

    std::vector<index_type> colsorted(std::distance(colbegin, colend));
    // shuffle the columns to get a sorted array of columns
    for (size_t i = 0;i < colsorted.size(); ++i) {
      colsorted[i] = *(colbegin + indices[i]);
    }
    // create the column indices
    std::unique_copy(colsorted.begin(), colsorted.end(),
                     std::inserter(jc, jc.end()));

    // resize the remaining indices
    cp.resize(jc.size() + 1);
    ir.resize(numvals);
    value.resize(numvals);


    // ok. cp is just a prefix sum of the colsorted array
    size_t cpidx = 0;
    for (size_t i = 1;i < colsorted.size(); ++i) {
      if (colsorted[i - 1] != colsorted[i]) {
        cp[cpidx + 1] = cp[cpidx];
        ++cpidx;
      }
      ++cp[cpidx];
    }
    cp[cp.size() - 1] = numvals;
    // ir and value are just shuffles of row and val
    for (size_t i = 0;i < ir.size(); ++i) {
      ir[i] = *(rowbegin + indices[i]);
      value[i] = *(valbegin + indices[i]);
    }
  }

  /**
   * Returns an iterator to the beginning of a list of columns
   */
  typename std::vector<index_type>::const_iterator col_list_begin() const {
    return jc.begin();
  }

  /**
   * Returns an iterator to the end of a list of columns
   */
  typename std::vector<index_type>::const_iterator col_list_end() const {
    return jc.end();
  }

  /**
   * Returns a column pointed to by a column iterator in the range
   * col_list_begin() to col_list_end();
   */
  entry_list get_column(typename std::vector<index_type>::const_iterator jciter) {
    if (jciter != jc.end()) {
      size_t c = std::distance(jc.begin(), jciter);
      size_t len = cp[c + 1] - cp[c];
      return entry_list((*jciter),
                        value.begin() + cp[c],
                        ir.begin() + cp[c],
                        len);
    } else {
      return entry_list((index_type)(-1),
                        value.end(),
                        ir.end(),
                        0);
    }
  }

/**
   * Returns a particular column. A binary search is needed to
   * find the specific column.
   */
  entry_list get_column(index_type col) {
    typename std::vector<index_type>::iterator jciter =
        std::lower_bound(jc.begin(), jc.end(), col);
    if (jciter != jc.end() && (*jciter) == col) {
      return get_column(jciter);
    } else {
      return entry_list((index_type)(-1),
                        value.end(),
                        ir.end(),
                        0);
    }
  }

  value_type& data_at_offset(size_t i) {
    :qa
  }
};

template <typename IndexType, typename ValueType>
std::ostream& operator<<(std::ostream& out,
                      const dcsc_store<IndexType, ValueType>& store) {
  store.print(out);
  return out;
}

} // namespace graphlab

#endif
