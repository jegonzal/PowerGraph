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
  template<typename valuetype, uint32_t blocksize=1024>
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
       assign(first, last);
     }

     /// Destructor
     ~block_linked_list() { clear(); }

     template<typename InputIterator>
     void assign(InputIterator first, InputIterator last)  {
       if (head != NULL)
         clear();
       InputIterator iter = first;
       head = tail = NULL;
       blocktype* next = NULL;
       float id = 0;
       while (iter != last) {
         InputIterator end =  std::min(iter+blocksize, last);
         next = new blocktype(id);
         next->assign(iter, end);
         iter = end;
         if (head == NULL) { head = next; }
         if (tail != NULL) { 
           tail->_next = next;
         }
         tail = next;
         next = tail->_next;
         id = id + 1;
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
    * Insert value into the location of iter. Split the block when  
    * necessary.
    */
   void insert(iterator iter, const valuetype& val) {
     blocktype* blockptr = iter.get_blockptr(); 
     ASSERT_TRUE(blockptr != NULL);
     uint32_t offset = iter.get_offset(); 
     if (blockptr->is_full()) {
       blockptr->split();
     }
     blockptr->insert(val, offset);
     ++_size;
   }

   blocktype* add_block() {
     float id = (tail == NULL) ? 0 : (tail->id() +1);
     blocktype* cur = new blocktype(id);
     if (tail == NULL) {
       head = tail = cur;
     } else {
       tail->_next = cur;
       tail = cur;
     }
     return cur;
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

   //////////////////// Pretty print API ///////////////////////// 
   void print(std::ostream& out) {
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
       head = other.head;
       tail = other.tail;
       _size = other._size;
       other.head = other.tail = NULL;
       other._size = 0;
     }

     void clear() {
       blocktype* tmp;
       while(head != tail) {
         tmp = head;
         head = head->_next;
         delete tmp;
       }
       delete tail;
       head = tail = NULL;
       _size = 0;
     }

     void load(iarchive& iarc) {
       // TODO
     }

     void save(oarchive& oarc) const {
       // TODO
     }
   
   //////////////////// Private Data Member ///////////////////////// 
   private:
     blocktype* head;
     blocktype* tail;
     size_t _size;
  };
} // end of namespace
#endif
