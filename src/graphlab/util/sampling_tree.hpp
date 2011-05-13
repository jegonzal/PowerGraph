/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_SAMPLING_TREE_HPP
#define GRAPHLAB_SAMPLING_TREE_HPP

#include <vector>


#include <boost/integer.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>

#include <graphlab/util/synchronized_queue.hpp>


namespace graphlab {
  /**
   \ingroup util_internal
  Fast multinomial sampler
  */
  template <typename Data, typename Priority>
  class sampling_tree{
  private:

    typedef boost::lagged_fibonacci607 random_source;
    //typedef boost::minstd_rand0 random_source;
    typedef boost::uniform_01<random_source> distribution_type;
    std::vector<random_source> rngs;
    std::vector<distribution_type> distributions;
  
  
    // each element in the tree is the sum of the elements in the children
    std::vector<Priority> tree;  

  
    // here we store all the leaves.
    std::vector< synchronized_queue<std::pair<Data, Priority> > > leaves;
    std::vector<spinlock> locks;
    atomic<size_t> num_elements;

    size_t max_key_value;
  
    size_t tree_interior_size;
  
    Priority minpriority;
  
    size_t nextpow2(size_t val) {
      --val;
      val = val | (val >> 1);
      val = val | (val >> 2);
      val = val | (val >> 4);
      val = val | (val >> 8);
      val = val | (val >> 16);
#ifdef _LP64   
      val = val | (val >> 32);
#endif
      return val + 1; 
    }
  
    /// Returns the index of a leaf sampled proportionate to its priority.
    /// Returns false on failure
    bool sample_tree(size_t &leaf, size_t cpuid) {
      size_t loc = 0;
      while (loc < tree_interior_size) {
        // get the left and right priorities
        Priority left = tree[loc * 2 + 1];
        Priority right = tree[loc * 2 + 2];
        // if both are zero, the sample has failed. Return
        if (left + right == 0) return false;
        else if (right == 0) {
          // if right is zero, move left
          loc = loc * 2 + 1;
        } else if (left == 0) {
          // if left is zero move right
          loc = loc * 2 + 2;
        } else {
          // pick from a bernoulli trial
          Priority childsum = left + right;

          double rndnumber = distributions[cpuid]() ;
          // double rndnumber = rand();
          Priority sel = rndnumber * childsum;
          loc = loc * 2 + 1 + (sel > left);
        }
      }
      leaf = loc - tree_interior_size;
      return true;
    }
  
    /// Pops an entry on a particular leaf. Returns false on failure
    bool pop_leaf(size_t leaf, std::pair<Data, Priority> &entry) {
      locks[leaf].lock();
      if (leaves[leaf].safepop(&entry)) {
      
        if (leaves[leaf].size() > 0) {
          tree[leaf + tree_interior_size] -= entry.second;
          if (tree[leaf + tree_interior_size] < 0) 
            tree[leaf + tree_interior_size] = 0;
        } else {
          tree[leaf + tree_interior_size] = 0;
        }
      
        locks[leaf].unlock();
        return true;
      }
      locks[leaf].unlock();
      return false;
    }

    /// Propagates a cumulative sum update up the tree.
    void propagate_up(size_t leaf) {
      leaf += tree_interior_size;
      while (leaf != 0) {
        volatile Priority *sibling1 = &(tree[leaf]);
        // the binary here is equivalent to (+1) if leaf is odd and (-1) if 
        // leaf is even
        volatile Priority *sibling2 = &(tree[leaf + (leaf & 1)*2 - 1]);
      
        leaf = (leaf - 1) >> 1;
        volatile Priority *parent = &(tree[leaf]);
        // write to the parent. Use a concurrent write mechanism
        do {
          (*parent) = (*sibling1) + (*sibling2);     
        } while((*parent) != (*sibling1) + (*sibling2));
      }
    }

  public:

    sampling_tree(size_t maxkeyval,
                  Priority minimumpriority = 1E-10,
                  size_t ncpus = 16) : rngs(ncpus){
      for(size_t i = 0; i < rngs.size(); ++i) {
        rngs[i].seed(rand());
        distributions.push_back(distribution_type(rngs[i]));
      
      }
    
      max_key_value = maxkeyval;
      minpriority = minimumpriority;
      // create the leaves
      //leaves.resize(maxkeyval, synchronized_queue<std::pair<Data, Priority> >(1));
      leaves.resize(maxkeyval);
      locks.resize(maxkeyval);
      // create the tree. (full binary tree
      tree_interior_size = nextpow2(maxkeyval) - 1;
      size_t treesize = 2 * nextpow2(maxkeyval) - 1;
      tree.resize(treesize, Priority(0));
      num_elements.value = 0;
    }

    ~sampling_tree(){}


    size_t size() {
      return num_elements.value;
    }
  
    bool pop(size_t &key , Data &ret, Priority &pr, size_t cpuid = 0) {
      if (size() == 0) return false;
      std::pair<Data, Priority> dppair;
      while (true) {
        if (sample_tree(key, cpuid)) {
          if (pop_leaf(key, dppair)) {
            propagate_up(key);
            break;
          }
        }
        if (size() == 0) {
          return false;
        }
      }
      ret = dppair.first;
      pr = dppair.second;
      num_elements.dec();
      return true;
    }
  
    bool sample(size_t &k, size_t cpuid = 0) {
      if (size() == 0) return false;
      while (true) {
        if (sample_tree(k, cpuid)) {
          return true;
        } else if (size() == 0) {
          return false;
        }
      }
    }

    void push(const size_t& key, const Data& item, const Priority &pr) {
      num_elements.inc();
    
      locks[key].lock();
      leaves[key].push(std::make_pair(item, std::max(pr, minpriority)));
      tree[key + tree_interior_size] += pr;
      locks[key].unlock();  
    
      propagate_up(key);
    }
 
    void increase_priority(const size_t& key, const Priority &pr) {
      // try to add my priority into the leaf. use CAS to make sure the
      // add gets in
      locks[key].lock();
      tree[key + tree_interior_size] += pr;
      locks[key].unlock();  

      propagate_up(key);    
    }
  
    void print_tree() {
      for (size_t i = 0;i < std::min(tree.size(),size_t(1000)); ++i) {
        std::cout << tree[i] << " ";
      }
    }
  };

}
#endif

