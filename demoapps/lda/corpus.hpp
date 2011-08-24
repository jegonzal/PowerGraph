/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef CORPUS_HPP
#define CORPUS_HPP

#include <stdint.h>
#include <cassert>
#include <string>
#include <vector>


typedef uint32_t word_id_type;
typedef uint32_t doc_id_type;
typedef uint16_t topic_id_type;
typedef size_t   count_type;

#define NULL_WORD word_id_type(-1)
#define NULL_TOPIC topic_id_type(-1)


struct corpus {


  struct token {
    word_id_type word;
    doc_id_type doc;
    token(const word_id_type& word = 0, const doc_id_type& doc = 0) : 
      word(word), doc(doc) { }
  };
  size_t nwords;
  size_t ndocs;
  size_t ntokens;
  std::vector< token > tokens;
  std::vector<std::string> dictionary;
  std::vector< word_id_type > ntokens_in_doc;
  corpus(const std::string& dictionary_fname, 
              const std::string& counts_fname );
  void load_dictionary(const std::string& fname);
  void load_counts(const std::string& fname);

  void shuffle_tokens();
}; // end of corpus

std::ostream& operator<<(std::ostream& out, const corpus::token& tok);

/**
 * Randomly split the corpus into two separate sets
 */
void split(const corpus& base, corpus& c1, corpus& c2, double prop_c1);


#endif
