#ifndef GRAPHLAB_DCSC_STORE_ITERATORS_HPP
#define GRAPHLAB_DCSC_STORE_ITERATORS_HPP
#include <vector>
#include <boost/version.hpp>
#include <boost/iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
namespace graphlab {

namespace dcsc_store_impl {
  // A class of edge information. Used as value type of the entry_list.
  template <typename EntryIndexType, typename EntryValueType>
  class entry_type {
   private:
    typedef EntryIndexType entry_index_type;
    typedef EntryValueType entry_value_type;
    entry_index_type col;
    typename std::vector<entry_value_type>::iterator valiter;
    typename std::vector<entry_index_type>::const_iterator rowiter;

   public:
     // Cosntructors
     /** \brief Creates an empty iterator. */
     entry_type () { }

     /** \brief Creates an iterator at a specific edge. */
     entry_type(entry_index_type col,
               typename std::vector<entry_value_type>::iterator valiter,
               typename std::vector<entry_index_type>::const_iterator rowiter)
         : col(col), valiter(valiter), rowiter(rowiter) { }

    /** \brief Returns the source vertex id of the edge. */
    inline const entry_index_type& column() const {
      return col;
    }
    /** \brief Returns the target vertex id of the edge. */
    inline const entry_index_type& row() const {
      return (*rowiter);
    }
    inline const entry_value_type& value() const {
      return (*valiter);
    }
    inline entry_value_type& value() {
      return (*valiter);
    }
  };


  template <typename EntryIndexType, typename EntryValueType>
  class entry_iterator  {
   public:
     typedef EntryIndexType entry_index_type;
     typedef EntryValueType entry_value_type;
     typedef entry_type<entry_index_type, entry_value_type> real_entry_type ;

     typedef std::random_access_iterator_tag iterator_category;
     typedef real_entry_type    value_type;
     typedef ssize_t      difference_type;
     typedef real_entry_type*   pointer;
     typedef real_entry_type   reference;

   private:
    entry_index_type col;
    typename std::vector<entry_value_type>::iterator valiter;
    typename std::vector<entry_index_type>::const_iterator rowiter;

   public:
     // Cosntructors
     /** \brief Creates an empty iterator. */
     entry_iterator () { }

     /** \brief Creates an iterator at a specific edge. */
     entry_iterator (entry_index_type col,
                    typename std::vector<entry_value_type>::iterator valiter,
                    typename std::vector<entry_index_type>::const_iterator rowiter)
         : col(col), valiter(valiter), rowiter(rowiter) { }

     inline real_entry_type operator*() const  {
        return make_value();

     }
#if BOOST_VERSION < 105000
     typedef boost::detail::
         operator_arrow_result<real_entry_type, real_entry_type, real_entry_type*> arrow_type;
     inline typename arrow_type::type operator->() {
       return arrow_type::make(make_value());
     }
#else
     typedef typename boost::detail::
         operator_arrow_dispatch<real_entry_type, real_entry_type*>::result_type arrow_type;
     inline arrow_type operator->() {
       return arrow_type(make_value());
     }
#endif

     /** \brief Returns if two iterators point to the same edge. */
     inline bool operator==(const entry_iterator& it) const {
       return col == it.col &&
           rowiter == it.rowiter;
     }

     /** \brief Returns if two iterators don't point to the same edge. */
     inline bool operator!=(const entry_iterator& it) const {
       return !(*this == it);
     }

     /** \brief Increases the iterator. */
     inline entry_iterator& operator++() {
       ++valiter;
       ++rowiter;
       return *this;
     }

     /** \brief Increases the iterator. */
     inline entry_iterator operator++(int) {
       const entry_iterator copy(*this);
       operator++();
       return copy;
     }

     /** \brief Computes the difference of two iterators. */
     inline ssize_t operator-(const entry_iterator& it) const {
       return rowiter - it.rowiter;
     }

     /** \brief Returns a new iterator whose value is increased by i difference units. */
     inline entry_iterator operator+(difference_type i) const {
       return entry_iterator(col, valiter + i, rowiter + i);
     }

     /** \brief Increases the iterator by i difference units. */
     inline entry_iterator& operator+=(difference_type i) {
       valiter += i;
       rowiter += i;
       return *this;
     }

     /** \brief Generate the return value of the iterator. */
     inline real_entry_type make_value() {
       return real_entry_type(col, valiter, rowiter);
     }

     /** \brief Generate the return value of the iterator. */
     inline real_entry_type make_value() const {
       return real_entry_type(col, valiter, rowiter);
     }

  };






  template <typename EntryIndexType, typename EntryValueType>
  struct entry_list {
    typedef EntryIndexType entry_index_type;
    typedef EntryValueType entry_value_type;
    typedef entry_type<entry_index_type, entry_value_type> real_entry_type ;

    entry_index_type col;
    typename std::vector<entry_value_type>::iterator valbegin;
    typename std::vector<entry_value_type>::iterator valend;
    typename std::vector<entry_index_type>::const_iterator rowbegin;
    typename std::vector<entry_index_type>::const_iterator rowend;

    typedef entry_iterator<entry_index_type, entry_value_type> iterator;
    typedef entry_iterator<entry_index_type, entry_value_type> const_iterator;
    typedef real_entry_type value_type;

    entry_list() { }

    entry_list(entry_index_type col,
              typename std::vector<entry_value_type>::iterator valbegin,
              typename std::vector<entry_index_type>::const_iterator rowbegin,
              size_t count)
        : col(col), valbegin(valbegin), valend(valbegin + count),
        rowbegin(rowbegin), rowend(rowbegin + count) { }

    iterator begin() {
      return iterator(col, valbegin, rowbegin);
    }
    iterator end() {
      return iterator(col, valend, rowend);
    }

    const_iterator begin() const {
      return const_iterator(col, valbegin, rowbegin);
    }
    const_iterator end() const {
      return const_iterator(col, valend, rowend);
    }
  };



} // namespace dcsc_store_impl
} // namespace graphlab

#endif
