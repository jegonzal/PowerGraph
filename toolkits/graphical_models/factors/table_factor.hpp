/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef TABLE_FACTOR_HPP
#define TABLE_FACTOR_HPP

/**
 * This file defines an opaque interface (or facade) for the various table types.
 *
 *  \author Scott Richardson     10/2012
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cassert>
#include <cmath>

#include <iostream>
#include <typeinfo>
#include <algorithm>
#include <limits>
#include <vector>
#include <set>

#include <graphlab/serialization/serialization_includes.hpp>
#include "dense_table.hpp"
#include "sparse_table.hpp"



namespace graphlab {



/**
 * An abstraction implemented to manage the creation and deletion 
 * of the sparse/dense tables; so that RAII
 */
 // NOTE im not 100% this class is necessary, except to provide 
 // convinient access to certian operators like operator*(). if i
 // made table_base concrete, i might be able to use it wherever
 // i currently use this class. (it's not like this class is 
 // providing conversions between dense and sparse tables or anything.)
template<size_t MAX_DIM>
class table_factor {
  typedef dense_table<MAX_DIM>    dense_table_t;
  typedef sparse_table<MAX_DIM>   sparse_table_t;
  typedef table_base<MAX_DIM>     table_base_t;

public:
  //static const size_t MAX_DIM_ = MAX_DIM; 
 

  // REVIEW accessing this requires constructing a table_factor
  enum table_storage_t { nil, DENSE_TABLE, SPARSE_TABLE };

  table_factor() : _table_storage(table_factor::nil), _table(NULL) { }

  table_factor(table_storage_t storage) : 
      _table_storage(storage), _table(NULL) 
  { 
    alloc();
  }

  table_factor(const table_factor& other) :
      _table_storage(other._table_storage), _table(NULL) 
  {
    if(other._table == NULL) { // integrity check
      DCHECK_EQ(other._table_storage, table_factor::nil);
    } else {
      alloc();

      _table->deep_copy(*(other._table));
    }
  }

  table_factor(table_storage_t storage, typename table_base_t::const_ptr base) : 
      _table_storage(storage), _table(NULL) 
  {
    DCHECK_NE(base, NULL);
    DCHECK_NE(_table_storage, table_factor::nil);

    alloc();
    _table->deep_copy(*base);
  }

  table_factor(table_storage_t storage, const table_base_t& base) : 
      _table_storage(storage), _table(NULL) 
  {
    DCHECK_NE(_table_storage, table_factor::nil);

    alloc();
    _table->deep_copy(base);
  }

  table_factor(typename table_base_t::const_ptr base) :
      _table_storage(table_factor::nil), _table(NULL) 
  {  
    DCHECK_NE(base, NULL);

    determine_storage_t(*base);
    alloc();

    _table->deep_copy(*base);
  }

  table_factor(const table_base_t& base) : 
      _table_storage(table_factor::nil), _table(NULL) {

    determine_storage_t(base);
    alloc();

    _table->deep_copy(base);
  }

  ~table_factor() {
    delete _table;
  }

  // TODO provide iterators

public:
  table_factor& operator=(const table_factor& other) {
    if(this == &other) return *this;

    table_base_t* base = NULL;
    if(other._table == NULL) { // integrity check
      DCHECK_EQ(other._table_storage, table_factor::nil);
    } else {
      try {
        alloc_table(other._table_storage, &base);
        base->deep_copy(*(other._table));
      } catch ( ... ) {
        delete base;
        throw;
      }
    }
    _table_storage = other._table_storage;
    delete _table;
    _table = base;

    return *this;
  }
/*
  table_factor& operator+=(const table_factor& other) {
    DCHECK_NE(_table, NULL);
    DCHECK_NE(other._table, NULL);

    _table->plus_equals(*(other._table));
    return *this;
  }
*/
  table_factor& operator*=(const table_factor& other) {
    DCHECK_NE(_table, NULL);
    DCHECK_NE(other._table, NULL);

    _table->times_equals(*(other._table));
    return *this;
  }

  table_factor& operator/=(const table_factor& other) {
    DCHECK_NE(_table, NULL);
    DCHECK_NE(other._table, NULL);

    _table->divide_equals(*(other._table));
    return *this;
  }

  table_factor operator*(const table_factor& other) const {
    DCHECK_NE(_table, NULL);
    DCHECK_NE(other._table, NULL);
    
    // deep copy
    table_factor out(_table_storage);
    _table->times(*(other._table), *(out._table));
    
    return out;
  }

  table_factor operator/(const table_factor& other) const {
    DCHECK_NE(_table, NULL);
    DCHECK_NE(other._table, NULL);
    
    // deep copy
    table_factor out(_table_storage);
    _table->divide(*(other._table), *(out._table));
    
    return out;
  }

  table_base_t const * table() const { 
    DCHECK_NE(_table, NULL);

    return _table; 
  }

  table_base_t* table() { 
    DCHECK_NE(_table, NULL);

    return _table; 
  }

  void load(graphlab::iarchive& arc) {
    arc >> _table_storage; 
    alloc();
    if(_table_storage != table_factor::nil) arc >> *_table;
  }

  void save(graphlab::oarchive& arc) const {
    arc << _table_storage;
    if(_table_storage != table_factor::nil) arc << *_table;
  }

private:
  // NOTE take ownership of base
  table_factor(table_storage_t storage, table_base_t *const *const base) : 
      _table_storage(storage) {
    if(base != NULL) _table = *base;
    DCHECK_NE(_table, NULL);
  }

  void determine_storage_t(const table_base_t& base) {
    if( typeid(base) == typeid(dense_table_t) ) 
      _table_storage = table_factor::DENSE_TABLE;
    else if( typeid(base) == typeid(sparse_table_t) ) 
      _table_storage = table_factor::SPARSE_TABLE;
    else {
      _table_storage = table_factor::nil;
      // REVIEW should probably raise an exception
      std::cout << "ERROR: unknown table storage type. " << std::endl;
      ASSERT_TRUE(false);
    }
  }

  void alloc() {
    alloc_table(_table_storage, &(this->_table));
  }

  // REVIEW could probably use virtual constructors to avoid some of this
  void alloc_table(table_storage_t storage, table_base_t** base) {
    // try to avoid memory leaks
    DCHECK_EQ(*base, NULL);

    switch(storage) {
      case table_factor::DENSE_TABLE:
        *base = new dense_table_t();
        break;
      case table_factor::SPARSE_TABLE:
        *base = new sparse_table_t();
        break;
      case table_factor::nil:
      default:
        // i allow *this to be constructed with a table_factor::nil table storage 
        // type, so this path is possible 
        break;
    }
  }

public:
  table_storage_t table_storage() const { return _table_storage; } 

  friend std::ostream& operator<<(std::ostream& out,
                            const table_factor<MAX_DIM>& factor) {
    out << "Table Factor(" << factor._table_storage << "): " << *(factor._table);
    return out;
  }

private: 
  table_storage_t _table_storage; 
  table_base_t* _table;

}; // end of table_factor



} // end of namespace graphlab 

#endif // TABLE_FACTOR_HPP
