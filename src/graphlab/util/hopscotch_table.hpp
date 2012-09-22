#ifndef GRAPHLAB_UTIL_HOPSCOTCH_TABLE_HPP
#define GRAPHLAB_UTIL_HOPSCOTCH_TABLE_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <iterator>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>

#if __cplusplus < 201103L 
#include <boost/functional/hash.hpp>
#define _HOPSCOTCH_TABLE_DEFAULT_HASH boost::hash<T>

#else
#define _HOPSCOTCH_TABLE_DEFAULT_HASH std::hash<T>
#endif

namespace graphlab {


/**
  * This defines a hash table where each entry stores a
  * fixed data type T. The data type T should be <b>small</b>
  * and should preferably fit in a couple of words.
  * This hash table is not resizeable. Use the hopscotch_map
  * For a more general purpose table.
  *
  * Safe access is guaranteed if you restrict to the functions suffixed with
  * _sync.
  *
  * \tparam T The data type stored in the hash table
  * \tparam Synchronized Defaults to True. If True, locking is used to ensure
  *                      safe reads and writes to the hash table.
  *                      Even under "Synchronized", the only operations
  *                      which are safe for parallel access are all functions 
  *                      suffixed with "sync".
  * \tparam Hash The hash functor type. Defaults to std::hash<T> if C++11 is
  *              available. Otherwise defaults to boost::hash<T>
  * \tparam KeyEqual The functor used to identify object equality. Defaults to
  *                  std::equal_to<T>
  */
template <typename T, 
         bool Synchronized = true,
         typename Hash = _HOPSCOTCH_TABLE_DEFAULT_HASH, 
         typename KeyEqual = std::equal_to<T> >
class hopscotch_table {
  public:
    /// The data type stored in the table
    typedef T              value_type;
    typedef size_t                                   size_type;
    /// The type of the hasher object
    typedef Hash                                     hasher;
    /// The type of the equality tester
    typedef KeyEqual equality_function;
    /// A pointer to the data type stored in the table
    typedef value_type* pointer;
    /// A reference to the data type stored in the table
    typedef value_type& reference;
    /// A constant pointer to the data type stored in the table
    typedef const value_type* const_pointer;
    /// A constant reference to the data type stored in the table
    typedef const value_type& const_reference;

  private:
    /// The actual contents of the hash table
    struct element {
      bool hasdata: 1;      /// Whether this entry has data.
      uint32_t field: 31;   /// The hopscotch bitfield. Only 31 bits are useable
      T elem;  /// User data
      element():hasdata(false), field(0) { }
    } __attribute__((__packed__));
    
    std::vector<element> data; 
    std::vector<simple_spinlock> locks;

    hasher hashfun;
    equality_function equalfun;
    atomic<size_t> numel;
    size_t mask;

    /// Returns the next power of 2 of a value
    static uint64_t next_powerof2(uint64_t val) {
      --val;
      val = val | (val >> 1);
      val = val | (val >> 2);
      val = val | (val >> 4);
      val = val | (val >> 8);
      val = val | (val >> 16);
      val = val | (val >> 32);
      return val + 1;
    }

    /** Computes the hash of the data. And perturbs it
      * using either CRC32 or Jenkin's 32-bit mix
      */
    size_t compute_hash(const value_type& d) const {
      size_t state = hashfun(d);
#ifdef HAS_BUILTIN_CRC32
      return __builtin_ia32_crc32di(0, state); 
#else
    /*
     * Bob Jenkin's 32 bit integer mix function from
     * http://home.comcast.net/~bretm/hash/3.html
     */
      state += (state << 12);
      state ^= (state >> 22);
      state += (state << 4);
      state ^= (state >> 9);
      state += (state << 10);
      state ^= (state >> 2);
      state += (state << 7);
      state ^= (state >> 12);
      return state;
#endif
    }


    /// Returns the lock ID associated with a given array index
    static size_t associated_lock_id(size_t idx) {
      return idx / 32;
    }

    

  public:
    /**
     * Constructs a hopscotch table of a given length.
     *
     * \param len This rounded up to the next power of 2 will be used as 
     *            the length of the table. This table is not resizeable.
     * \param hashfun The hasher functor. Defaults to Hash()
     * \param equalfun A functor used to test for equality. Defaults to KeyEqual()
     */
    hopscotch_table(size_t len,
                    Hash hashfun = Hash(),
                    KeyEqual equalfun = KeyEqual()):
                                              data(next_powerof2(len) + 32), 
                                              locks(Synchronized ? data.size() / 32 + 1 : 0),
                                              hashfun(hashfun),
                                              equalfun(equalfun),
                                              numel(0),
                                              mask(data.size() - 32 - 1) {
    }

    /// Returns the hash function used by the hash table
    hasher hash_function() const {
      return hashfun;
    }

    /// Returns the equality function used by the hash table
    equality_function key_eq() const {
      return equalfun;
    }

    /**
      * A const iterator which allows iteration over the hash table
      * entries. Insertions may disrupt the iterator order. Deletions
      * invalidate the iterator.
      */
    struct const_iterator {
      typedef std::forward_iterator_tag iterator_category;
      typedef const typename hopscotch_table::value_type value_type;
      typedef size_t difference_type;
      typedef value_type* pointer;
      typedef value_type& reference;

      friend class hopscotch_table;

      const hopscotch_table* ptr;
      typename std::vector<element>::const_iterator iter;

      const_iterator():ptr(NULL) {}

      const_iterator operator++() {
        ++iter;
        while(iter != ptr->data.end() && !iter->hasdata) {
          ++iter;
        } 
        return *this;
      }

      const_iterator operator++(int) {
        iterator cur = *this;
        ++(*this);
        return cur;
      }


      reference operator*() {
        return iter->elem; 
      }

      pointer operator->() {
        return &(iter->elem); 
      }

      bool operator==(const const_iterator it) const {
        return ptr == it.ptr && iter == it.iter;
      }

      bool operator!=(const const_iterator iter) const {
        return !((*this) == iter);
      }


    private:
      const_iterator(const hopscotch_table* table, 
          typename std::vector<hopscotch_table::element>::const_iterator iter): 
        ptr(table), iter(iter) { }
    };



    /**
      * A const iterator which allows iteration over the hash table
      * entries. Insertions may disrupt the iterator order. Deletions
      * invalidate the iterator.
      */
    struct iterator {
      typedef std::forward_iterator_tag iterator_category;
      typedef typename hopscotch_table::value_type value_type;
      typedef size_t difference_type;
      typedef value_type* pointer;
      typedef value_type& reference;

      friend class hopscotch_table;

      hopscotch_table* ptr;
      typename std::vector<element>::iterator iter;

      iterator():ptr(NULL) {}


      operator const_iterator() const {
        const_iterator it(ptr, iter);
        return it;
      }

      iterator operator++() {
        ++iter;
        while(iter != ptr->data.end() && !iter->hasdata) {
          ++iter;
        } 
        return *this;
      }

      iterator operator++(int) {
        iterator cur = *this;
        ++(*this);
        return cur;
      }


      reference operator*() {
        return iter->elem; 
      }

      pointer operator->() {
        return &(iter->elem); 
      }

      bool operator==(const iterator it) const {
        return ptr == it.ptr && iter == it.iter;
      }

      bool operator!=(const iterator iter) const {
        return !((*this) == iter);
      }


    private:
      iterator(hopscotch_table* table, 
          typename std::vector<element>::iterator iter):
        ptr(table), iter(iter) { }
    };

    /** Standard insert iterator. Writing into this iterator
     *  will cause insertions to occur. It is however, recommended
     *  that the insert() operation be used instead of the insert_iterator
     *  since the insert_iterator silently fails on insert failure.
     */
    struct insert_iterator{
      hopscotch_table* cmap;
      typedef std::forward_iterator_tag iterator_category;
      typedef typename hopscotch_table::value_type value_type;

      insert_iterator(hopscotch_table* c):cmap(c) {}
      
      insert_iterator operator++() {
        return (*this);
      }
      insert_iterator operator++(int) {
        return (*this);
      }

      insert_iterator& operator*() {
        return *this;
      }
      insert_iterator& operator=(const insert_iterator& i) {
        cmap = i.cmap;
        return *this;
      }
      
      insert_iterator& operator=(const value_type& v) {
        cmap->insert(v);
        return *this;
      }
    };
   
  private: 
    /**
     *  Searches for a target entry and overwrites if it exists
     */
    template <bool IsSynchronized>
    iterator try_find_and_overwrite(const value_type& newdata,
                                    size_t target,
                                    bool overwrite) {
       // find the next empty entry 
      size_t lockid = associated_lock_id(target);
      
      if (IsSynchronized) {
        locks[lockid + 1].lock();
        locks[lockid].lock();
      }
      iterator iter = find_impl(newdata, target);
      if (iter != end() && overwrite) {
        iter.iter->elem = newdata;
      }
      if (IsSynchronized) {
        locks[lockid].unlock();
        locks[lockid + 1].unlock();
      }
      return iter;
    }

    /** Insert logic. If IsSynchronized is set, locks are used
     *  If overwrite is set, it will additionally check for existance
     *  of the entry and overwrite if it exists.
     * Iterator is not going to be necessarily valid under parallel access.
     */
    template <bool IsSynchronized>
    iterator insert_impl(const value_type& newdata, bool overwrite = true) {
      // find the next empty entry 
      size_t target = compute_hash(newdata) & mask;

      iterator ret = try_find_and_overwrite<IsSynchronized>(newdata, 
                                                            target, 
                                                            overwrite);
      if (ret != end()) return ret;

      // search for a place to stick it into
      bool found = false;
      size_t shift_target = target;
      // let max range is 31 * 20
      size_t limit = std::min(data.size(), target + 31 * 20);
      size_t lockid = 0;
      for (;shift_target < limit; shift_target++) {
        if (data[shift_target].hasdata == false) {
          lockid = associated_lock_id(shift_target);
          if (IsSynchronized) locks[lockid].lock();
          // double check
          if (data[shift_target].hasdata == false) {
            // yup still true.
            // we got an empty value.
            // quit the search
            found = true;
            break;
          } else {
            // nope. not empty anymore
            // unlock and continue
            if (IsSynchronized) locks[lockid].unlock();
          }
        }
      }

      if (!found) {
        // failed to find a place to put this value.
        return iterator(this, data.end());
      }

      // while the shift target is out of range
      while(shift_target - target >= 31) { 
        // search backwards
        // we would like to jump as far as possible
        // find an hash entry whose field placed something
        // between here and the shift target
        // and move it to the shift target.
        // for i = 31 to 1
        found = false;
        // lock one before the current lockid if available
        if (IsSynchronized && lockid > 0) locks[lockid - 1].lock();

        for (size_t i = 30; i >= 1; --i) {
          size_t r;
          if (data[shift_target - i].field) {
            r = __builtin_ctz(data[shift_target - i].field);
            if (r <= i) {
              // shift
              size_t new_shift_target = shift_target - i + r;
              assert(data[new_shift_target].hasdata);
              data[shift_target].elem = data[new_shift_target].elem;
              data[shift_target].hasdata = true;
              data[new_shift_target].hasdata = false;
              data[new_shift_target].elem = T();

              // unset the bit for r and set the bit for i
              data[shift_target - i].field = 
                (data[shift_target - i].field & ~((uint32_t)1 << r)) 
                 | ((uint32_t)1 << i);
              shift_target = new_shift_target;
              found = true;
              break;
            }
          }
        }

        if (!found) {
          // release all the locks acquired
          if (IsSynchronized) {
            locks[lockid].unlock();
            if (lockid > 0) locks[lockid - 1].unlock();
          }
          return iterator(this, data.end());
        }
        else {

          if (IsSynchronized) {
            // ok. depending on how far we went. we need to 
            // unlock one of lockid or lockid - 1
            size_t newlockid = associated_lock_id(shift_target);
            assert(newlockid == lockid || newlockid == lockid - 1);
            if (newlockid == lockid) {
              if (lockid > 0) locks[lockid - 1].unlock();
            }
            else if (newlockid == lockid - 1) {
              locks[lockid].unlock();
            }
            lockid = newlockid;
          }
        }
      }
      // insert and return
      // we need to lock ID - 1 so as to ensure intersection with the hash target
      if (IsSynchronized && lockid > 0) locks[lockid - 1].lock();
      data[shift_target].elem = newdata;
      data[target].field |= (1 << (shift_target - target));
      data[shift_target].hasdata = true;
      if (IsSynchronized) {
        ++numel;
        if (lockid > 0) locks[lockid - 1].unlock();
        locks[lockid].unlock();
      }
      else {
        ++numel.value;
      }
      return iterator(this, data.begin() + shift_target);
    }


    /**
      * Searches for an entry and returns an iterator to the entry.
      * The hash entry of the key is provided.
      * KeyEqual will be used to identify if an entry matches the request.
      * return end() on failure.
      */
    const_iterator find_impl(const value_type& key, size_t target) const {
      uint32_t field = data[target].field;
      while (field > 0) {
        int r = __builtin_ctz(field);
        if (data[target + r].hasdata && 
            key_eq()(data[target + r].elem, key)) {
          return const_iterator(this, data.begin() + target + r);
        }
        else {
          // mask out the current bit and try again.
          field &= ~(((uint32_t)1 << r));
        }
      }
      return const_iterator(this, data.end());
    }


    iterator find_impl(const value_type& key, size_t target) {
      const_iterator iter = ((const hopscotch_table*)(this))->find_impl(key, target);
      return iterator(this, data.begin() + (iter.iter - data.begin()));
    }

  public:
    /**
      * Inserts an entry into the array.
      * Returns an iterator to the just inserted data on success.
      * If the entry already exists, it will be overwritten.
      * Returns end() on failure.
      */
    iterator insert(const value_type& newdata) {
      return insert_impl<false>(newdata);
    }

    /**
      * Inserts an entry into the array.
      * Returns an iterator to the just inserted data on success.
      * This function check if the entry already exists, if it does, 
      * do nothing
      * Returns end() on failure.
      */
    iterator insert_do_not_overwrite(const value_type& newdata) {
      return insert_impl<false>(newdata, false);
    }



    /**
      * Searches for an entry and returns an iterator to the entry.
      * KeyEqual will be used to identify if an entry matches the request.
      * return end() on failure.
      */
    const_iterator find(const value_type& key) const {
      size_t target = compute_hash(key) & mask;
      return find_impl(key, target);
    }

    /**
      * Searches for an entry and returns an iterator to the entry.
      * KeyEqual will be used to identify if an entry matches the request.
      * return end() on failure.
      */
    iterator find(const value_type& key) {
      const_iterator iter = ((const hopscotch_table*)(this))->find(key);
      return iterator(this, data.begin() + (iter.iter - data.begin()));
    }

    


   /**
    * Erases an entry pointed to by an iterator.
    */
    bool erase(iterator iter)  {
      if (iter.iter == data.end()) return false;
      assert(iter.iter->hasdata);
      size_t target = compute_hash(iter.iter->elem) & mask;
      size_t offset = iter.iter - (data.begin() + target);
      assert(offset < 31);
      --numel.value;
      iter.iter->hasdata = false;
      iter.iter->elem = value_type();
      data[target].field &=  ~((uint32_t)1 << offset);
      return true;
    }

    /// Erases a entry matching a given value.
    bool erase(const value_type& key) {
      return erase(find(key));
    }

    /// Returns an iterator to the start of the table
    iterator begin() {
      // find the first which is not empty
      typename std::vector<element>::iterator iter = data.begin();      
      while (iter != data.end() && !iter->hasdata) {
        ++iter;
      }
      return iterator(this, iter);
    }

    /// Returns an iterator to the start of the table
    const_iterator begin() const {
      // find the first which is not empty
      typename std::vector<element>::iterator iter = data.begin();      
      while (iter != data.end() && !iter->hasdata) {
        ++iter;
      }
      return const_iterator(this, iter);
    }

    /// Returns an iterator to the end of the table
    iterator end() {
      return iterator(this, data.end());
    }

    /// Returns an iterator to the end of the table
    const_iterator end() const {
      return const_iterator(this, data.end());
    }

    /// Returns 1 if the table contains a given element. 0 otherwise.
    size_t count(const value_type& v) const {
      return find(v) != end();
    }

    /// Returns true if the table contains a given element. false otherwise.
    bool contains(const value_type& v) const {
      return find(v) != end();
    }

    /// Returns the number of elements in the table
    size_t size() const {
      return numel;
    }

    /// Returns the capacity of the table
    size_t capacity() const {
      return data.size();
    }

    float load_factor() const {
      return float(size()) / capacity();
    }

    // now for the safe accessors

    /// Returns the size of the hash table. Safe under parallel access.
    size_t size_sync() const {
      return numel;
    }

    /// Returns the capacity of the table
    size_t capacity_sync() const {
      return data.size();
    }

    float load_factor_sync() const {
      return float(size_sync()) / capacity_sync();
    }

    hopscotch_table& operator=(const hopscotch_table& other) {
      data = other.data;
      locks.resize(other.locks.size());
      hashfun = other.hashfun;
      equalfun = other.equalfun;
      numel = other.numel;
      mask = other.mask;
    }


    /** Inserts an element into the hash table. Safe under parallel access.
      * if t already exists, it will be overwritten
      */
    bool put_sync(const T& t) {
      // since data is not resizeable, 
      // data.end() is always valid.
      return insert_impl<Synchronized>(t).iter != data.end();
    }


    /** Inserts an element into the hash table. Safe under parallel access.
      * if t already exists, nothing will happen
      */
    bool put_do_not_overwrite_sync(const T& t) {
      // since data is not resizeable, 
      // data.end() is always valid.
      return insert_impl<Synchronized>(t, false).iter != data.end();
    }


    /** If the argument is found in the hash table,
     *  return {true, V} where V is the hash table content matching the argument.
     *  Otherwise {false, T()} is returned.
     *  KeyEqual() is used to compare entries.
     *  Safe under parallel access.
     */
    std::pair<bool, T> get_sync(const T& t) const {
      // fast path. Try to get it without locking
      const_iterator iter = find(t);
      if (iter != end()) {
        // take a snapshot of the data
        element e = *(iter.iter);
        if (e.hasdata && key_eq()(e.elem, t)) {
          return std::make_pair(true, e.elem);
        }
      }

      if (!Synchronized) {
        return std::make_pair(false, T());
      }

      // slow path. lock 
      size_t target = compute_hash(t) & mask;
      size_t lockid = associated_lock_id(target);
      locks[lockid + 1].lock();
      locks[lockid].lock();
      iter = find_impl(t, target);
      if (iter != end()) {
        element e = *(iter.iter);
        locks[lockid].unlock();
        locks[lockid + 1].unlock();
        assert(e.hasdata && key_eq()(e.elem, t));
        return std::make_pair(true, e.elem);
      }
      else {
        locks[lockid].unlock();
        locks[lockid + 1].unlock();
        return std::make_pair(false, T());
      }
    }

    /**
      * Returns true if the erasure was successful. 
      * If false, the entry was not found in the hash table
      */
    bool erase_sync(const T& t) {
      if (!Synchronized) {
        return erase(t);
      }

      // we need to find it first
      size_t target = compute_hash(t) & mask;
      size_t lockid = associated_lock_id(target);
      // acquire locks around the target
      locks[lockid + 1].lock();
      locks[lockid].lock();
      iterator iter = find(t);
      if (iter == end()) {
        locks[lockid].unlock(); 
        locks[lockid + 1].unlock();
        return false;
      }
      // now lets erase it
      assert(iter.iter->hasdata);
      size_t offset = iter.iter - (data.begin() + target);
      assert(offset < 31);
      --numel;
      iter.iter->hasdata = false;
      iter.iter->elem = value_type();
      data[target].field &=  ~((uint32_t)1 << offset);
      locks[lockid + 1].unlock();
      locks[lockid].unlock();
      return true;
    }
};

} // graphlab
#endif
