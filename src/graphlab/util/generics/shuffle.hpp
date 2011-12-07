#ifndef GRAPHLAB_INPLACE_SHUFFLE_HPP
#define GRAPHLAB_INPLACE_SHUFFLE_HPP
#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#include <omp.h>


namespace graphlab {
/**
 * Shuffles a random access container inplace such at
 * newcont[i] = cont[targets[i]]
 * targets must be the same size as the container
 * Both the container and the targets vector will be modified.
 */
template <typename Iterator, typename sizetype>
void inplace_shuffle(Iterator begin,
                     Iterator end, 
                     std::vector<sizetype> &targets) {
  size_t len = std::distance(begin, end);
  assert(len == targets.size());
  
  for (size_t i = 0;i < len; ++i) {
    // begin the permutation cycle
    if (i != targets[i]) {
      typename std::iterator_traits<Iterator>::value_type v = *(begin + i);
      size_t j = i;
      while(j != targets[j]) {
        size_t next = targets[j];
        if (next != i) {
          *(begin + j) = *(begin + next);
          targets[j] = j;
          j = next;
        }
        else {
          // end of cycle
          *(begin + j) = v;
          targets[j] = j;
          break;
        }
      }
    }
  }
}

/**
 * Shuffles a random access container inplace such at
 * newcont[i] = cont[targets[i]]
 * targets must be the same size as the container
 */
template <typename Container, typename sizetype>
void outofplace_shuffle(Container &c,
                        const std::vector<sizetype> &targets) {  
  Container result(targets.size());

#pragma omp parallel for
  for (ssize_t i = 0;i < ssize_t(targets.size()); ++i) {
    result[i] = c[targets[i]];
  }
  std::swap(c, result);
}

}
#endif
