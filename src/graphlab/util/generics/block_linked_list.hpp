#ifndef GRAPHLAB_BLOCK_LINK_LIST_HPP
#define GRAPHLAB_BLOCK_LINK_LIST_HPP
#include <graphlab/util/generics/dynamic_block.hpp>
#include <graphlab/logger/assertions.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <stdint.h>
#include <algorithm>

namespace graphlab {
  template<typename valuetype, uint32_t blocksize=64>
  class block_linked_list {
   public:
     typedef dynamic_block<valuetype, blocksize> blocktype;

   //////////////////// Constructors ///////////////////////// 
   public:
     /// Construct empty list
     block_linked_list() : head(NULL), tail(NULL), _size(0) { }

     /// Construct list from container 
     template<typename InputIterator>
     block_linked_list(InputIterator first, InputIterator last) { 
       if (first == last) 
         return;
       assign(first, last);
     }

     /// Destructor
     ~block_linked_list() { clear(); }

     template<typename InputIterator>
     void assign(InputIterator first, InputIterator last)  {
       if (_size > 0)
         clear();
       InputIterator iter = first;

       blocktype* current = head;
       if (current == NULL) {
         current = head = tail = new blocktype(0);
       }
       float id = 0;
       while (iter != last) {
         InputIterator end =  std::min(iter+blocksize, last);
         current->assign(iter, end);
         iter = end;
         id = id + 1;
         if (iter != last) {
           current = new blocktype(id);
           tail->_next = current;
           tail = current;
         }
       }
       _size = last - first;
     }

     /// Returns the size of the list
     size_t size() const { return _size; }

     static const size_t get_blocksize() {
       return blocksize;
     }

   //////////////////// Iterator API ///////////////////////// 
   public:
     template <typename value>
     class value_iterator :
         public boost::iterator_facade<value_iterator<value>, value,
                                       boost::random_access_traversal_tag> {
     private:
        struct enabler {};

     public:
        value_iterator(blocktype* blockptr, uint32_t offset) : 
            blockptr(blockptr), offset(offset) { }

        template <typename othervalue>
            value_iterator(value_iterator<othervalue> const& other,
                        typename boost::enable_if<
                        boost::is_convertible<othervalue*,value*>,
                        enabler>::type = enabler()) 
            : blockptr(other.blockptr), offset(other.offset) { }

     public: // returns block iterator
        blocktype*& get_blockptr() {
          return blockptr;
        }
        uint32_t& get_offset() {
          return offset;
        }

      private: // core access functions
        friend class boost::iterator_core_access;
        template <class> friend class value_iterator;

        void increment() { 
          if ((offset+1) < blockptr->size()) {
            ++offset;
          } else {
            blockptr = blockptr->_next;
            offset = 0;
          }
        }
        template <typename othervalue>
        bool equal(value_iterator<othervalue> const& other) const {
          return (blockptr == other.blockptr) && (offset == other.offset);
        }
        value& dereference() const { 
          return blockptr->values[offset];
        }
        void advance(int n) {
          size_t dist = n+offset;
          while(dist >= blockptr->size() && blockptr != NULL) {
            dist -= blockptr->size();
            blockptr = blockptr->next();
          } 
          if (blockptr == NULL) {
            offset = 0;
          } else {
            offset = dist;
          }
        } 
        ptrdiff_t distance_to(const value_iterator& other) const {
          ptrdiff_t dist = 0;
          if (blockptr == other.blockptr) {
            dist = (ptrdiff_t)other.offset - offset;
          } else {
            // determine the moving direction: forward if this is before other;
            bool move_forward = (other.blockptr == NULL) || (!(blockptr == NULL) && blockptr->id() < other.blockptr->id());
            if (move_forward) {
              blocktype* cur = blockptr;
              while (cur != other.blockptr) {
                dist += cur->size();
                cur = cur->next();
              }
              dist = dist + (ptrdiff_t)other.offset - offset;
            } else {
              // this after other
              blocktype* cur = other.blockptr;
              while (cur != blockptr) {
                dist += cur->size();
                cur = cur->next();
              }
              dist = -(dist + (ptrdiff_t)offset - other.offset);
            }
          }
          return dist;
        }
      private:
        blocktype* blockptr;
        uint32_t offset;
     }; // end of value iterator

     typedef value_iterator<valuetype> iterator; 
     typedef value_iterator<valuetype const> const_iterator; 

     iterator begin() {
       return iterator(head, 0);
     }

     iterator end() {
       return iterator(NULL, 0);
     }

     const_iterator begin() const {
       return const_iterator(head, 0);
     }

     const_iterator end() const {
       return const_iterator(NULL, 0);
     }

   //////////////////// Insertion API ///////////////////////// 
   /*
    * Insert value into the location of iter.
    * Split the block when necessary.
    * Returns iterator to the new value.
    */
   iterator insert(iterator iter, const valuetype& val) {
     iterator ins_iter = get_insert_iterator(iter);
     blocktype* ins_ptr = ins_iter.get_blockptr();
     uint32_t offset = ins_iter.get_offset();
     ASSERT_TRUE(ins_ptr != NULL);
     ++_size;
     if (ins_ptr->is_full()) {
       ins_ptr->split();
       if(ins_ptr == tail) {
         tail = ins_ptr->next();
       }
       if (offset >= blocksize/2) {
         ins_ptr = ins_ptr->next();
         offset -= (blocksize/2);
       }
     }
     ins_ptr->insert(val,offset);
     return iterator(ins_ptr, offset);
   }
   
   /**
    * Insert a range of values into the position of the given iterator.
    * Will create new blocks after the given block when necessary.
    *
    * \note 
    * This operation will NOT affect the blocks after the given block.
    *
    * Returns the begin and end iterator to the new elements. 
    *
    *
    * |x1,x2,y1,y2,y3,y4 , _ , _ , _ , _| -> | ... | -> |...|
    *        ^ 
    *        p
    * blocksize = 10
    * iterator: blockptr = p, offset = 2
    */
   template<typename InputIterator>
   std::pair<iterator, iterator> 
     insert(iterator iter, InputIterator first, InputIterator last) {

     typedef std::pair<iterator, iterator> ret_type;
     const size_t len = last - first;
     if (len == 0) return ret_type(iter, iter);

     iterator ins_iter = get_insert_iterator(iter);

     // Pointers to the block of the insertion point.
     blocktype* ibegin_ptr = ins_iter.get_blockptr(); 

     size_t nx = ins_iter.get_offset();
     size_t ny = ibegin_ptr->size()-nx;

     // save y 
     valuetype* swap = (valuetype*)malloc((ny)*sizeof(valuetype));
     memcpy(swap, &(ibegin_ptr->values[nx]), ny*sizeof(valuetype)) ;

     // remove y temporarily
     ibegin_ptr->_size -= ny;
     _size -= ny;

     // Insert new elements, keep begin and end iterators
     ret_type iter_pair = append_to_block(ibegin_ptr, first, last);

     // add y back 
     blocktype* iend_ptr = iter_pair.second.get_blockptr();
     ret_type iter_pair2 = append_to_block(iend_ptr, swap, swap+ny); 

     // Collect begin and end iterators
     iterator begin_ins_iter = iter_pair.first;
     iterator end_ins_iter = iter_pair2.first;

     if (end_ins_iter.get_offset() == 
         end_ins_iter.get_blockptr()->size()) {
       end_ins_iter.get_blockptr() = end_ins_iter.get_blockptr()->next();
       end_ins_iter.get_offset() = 0;
     }
     return ret_type(begin_ins_iter, end_ins_iter);
   }

   /**
    * Insert a range of values into the end of the given block.
    * Will create new blocks after the given block when necessary.
    *
    * \note 
    * This operation will NOT affect the blocks after the given block.
    *
    * Returns the begin and end iterator to the new elements. 
    */
  template<typename InputIterator>
  std::pair<iterator, iterator> 
     append_to_block(blocktype* ibegin_ptr, InputIterator first, InputIterator last) {
     ASSERT_TRUE(ibegin_ptr != NULL);

     const size_t len = last-first;
     blocktype* iend_ptr = NULL;
     uint32_t ibegin_offset, iend_offset;
     
     size_t nold = ibegin_ptr->size();
     size_t spaceleft = (blocksize - nold); 
     size_t nnew = std::min(len, spaceleft);

     // Fill in the rest of the block
     ASSERT_TRUE(nold+nnew <= blocksize);
     std::copy(first, first+nnew, &(ibegin_ptr->values[nold]));
     ibegin_ptr->_size += nnew;
     ibegin_offset = nold;

     iend_ptr = ibegin_ptr;
     iend_offset = nold+nnew;

     // creates a block chain for remaining elements
     if (len > spaceleft) {
       blocktype* current = insert_block(ibegin_ptr); 
       InputIterator iter = first + spaceleft; 
       while (iter != last) {
         InputIterator end =  std::min(iter+blocksize, last);
         current->assign(iter, end);
         iter = end;
         if (iter != last) {
           current = insert_block(current); 
         }
       }
       iend_ptr = current;
       iend_offset = iend_ptr->size(); 
     }

     _size += len;

     return std::pair<iterator,iterator>(
         iterator(ibegin_ptr, ibegin_offset)
         ,iterator(iend_ptr, iend_offset));
   }

   //////////////////// Block Access API ///////////////////////// 
   /*
    * Returns the nth block in the list. Linear time.
    */
   blocktype* get_block(size_t n) {
     size_t i = 0;
     blocktype* cur = head;
     while (cur != NULL && i < n) {
       cur = cur->_next;
       ++i;
     }
     return cur;
   }

   blocktype* insert_block(blocktype* ins) {
     float id = (ins==tail) ? 
         tail->id() +1 : (ins->id()+ins->next()->id())/2 ;
     blocktype* ret = new blocktype(id);
     ret->_next = ins->next();
     ins->_next = ret;
     if (ins == tail) {
       tail = ret;
     }
     return ret;
   }

   blocktype* append_block() {
     return insert_block(tail);
   }

   size_t num_blocks() const {
     if (head == NULL) 
       return 0;
     size_t ret = 1;
     blocktype* ptr = head;
     while (ptr != tail) {
       ptr = ptr->next();
       ++ret;
     }
     return ret;
   }

   //////////////////// Pretty print API ///////////////////////// 
   void print(std::ostream& out) const {
     blocktype* cur = head;
     while (cur != NULL) {
       cur->print(out);
       out << "-> ";
       cur = cur->_next;
     }
     out << "||" << std::endl;
   }
   
   //////////////////// Read Write ///////////////////////// 
     void swap (block_linked_list& other) {
       clear();
       delete head;
       head = other.head;
       tail = other.tail;
       _size = other._size;
       other.head = other.tail = NULL;
       other._size = 0;
     }

     void clear() {
       blocktype* tmp = NULL;
       while(head != tail) {
         ASSERT_TRUE(head != NULL);
         tmp = head;
         head = head->_next;
         delete tmp;
       }
       delete head;
       head = tail = NULL;
       _size = 0;
     }

     void load(iarchive& iarc) {
       // TODO
     }

     void save(oarchive& oarc) const {
       // TODO
     }
   

   //////////////////// Helper Function ///////////////////////// 
   iterator get_insert_iterator(iterator iter) {
     bool is_end = (iter == end());
     if (is_end) {
       if (tail == NULL) {
         head = tail = new blocktype(0);
       } else if (tail->is_full()) {
         append_block();
       } 
       iter.get_blockptr() = tail;
       iter.get_offset() = tail->size();
     }
     return iter;
   }

   //////////////////// Private Data Member ///////////////////////// 
   private:
     blocktype* head;
     blocktype* tail;
     size_t _size;
  };
} // end of namespace
#endif

     // } 
     // else {
     //   blocktype* splitblk = insert_block(ibegin_ptr);
     //   size_t nfirst = (nold+nnew)/2; 
     //   size_t nsecond = (nold+nnew-nfirst);
     //   ibegin_ptr->_size = nfirst;
     //   splitblk->_size = nsecond;
     //   if (nold < nfirst) {
     //     size_t padsize = nfirst - nold;
     //     std::copy(first, first+padsize, &(ibegin_ptr->values[nold]));
     //     std::copy(first+padsize, first+nnew, splitblk->values);
     //     ibegin_offset = nold;
     //   } else {
     //     size_t padsize = nold-nfirst;
     //     memcpy(splitblk->values, &(ibegin_ptr->values[nfirst]), padsize*sizeof(valuetype));
     //     std::copy(first, first+nnew, &(splitblk->values[padsize]));
     //     ibegin_ptr = splitblk;
     //     ibegin_offset = padsize;
     //  }
       // if (!add_new_blocks) {
       //    iend_ptr = splitblk;
       //    iend_offset = nsecond;
       // }
     // }

