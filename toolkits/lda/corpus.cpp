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


#include <fstream>
#include <algorithm>

#include <graphlab.hpp>

#include "corpus.hpp"



std::ostream& operator<<(std::ostream& out, const corpus::token& tok) {
  return out << "(" << tok.word << ", " << tok.doc << ")";
}


corpus::corpus(const std::string& dictionary_fname, 
               const std::string& counts_fname ) : 
  nwords(0), ndocs(0), ntokens(0), ntokens_in_doc(50000) {
  dictionary.reserve(20000);
  tokens.reserve(100000);
  load_dictionary(dictionary_fname);
  load_counts(counts_fname);
}



corpus::corpus(const std::string& dictionary_fname, 
               const std::string& counts_fname,
               const size_t procid, const size_t nprocs) : 
  nwords(0), ndocs(0), ntokens(0), ntokens_in_doc(50000) {
  dictionary.reserve(20000);
  tokens.reserve(100000);
  load_dictionary(dictionary_fname);
  load_counts(counts_fname, procid, nprocs);
}



void corpus::load_dictionary(const std::string& fname) {
  std::ifstream fin(fname.c_str());
  std::string str;
  while(fin.good()) {
    fin >> str;
    if(fin.good()) {
      dictionary.push_back(str);
      nwords++;
    }
  }
  fin.close();
}

void corpus::load_counts(const std::string& fname) {
  std::ifstream fin(fname.c_str());    
  while(fin.good()) {
    // Read a collection of tokens
    size_t word=0, doc=0, count=0;
    fin >> word >> doc >> count;
    if(fin.good()) {
      assert(word > 0 && doc > 0 && count > 0);
      // Decrement to c-style from matlab style
      word--; doc--;
      // update the doc counter
      ndocs = std::max(ndocs, doc + 1);
      // Assert valid word
      assert(word < nwords);
      // Update the words in document counter
      ntokens_in_doc[doc] += count;
      // Add all the tokens
      corpus::token tok; tok.word = word; tok.doc = doc;
      for(size_t i = 0; i < count; ++i) tokens.push_back(tok);
    }
  }
  fin.close();
  ntokens = tokens.size();
}


void corpus::load_counts(const std::string& fname,
                         const size_t procid,
                         const size_t nprocs) {
  std::ifstream fin(fname.c_str());    
  while(fin.good()) {
    // Read a collection of tokens
    size_t word=0, doc=0, count=0;
    fin >> word >> doc >> count;
    if(fin.good()) {
      assert(word > 0 && doc > 0 && count > 0);
      // Decrement to c-style from matlab style
      word--; doc--;
      if( (doc + procid) % nprocs == 0  ) {
        // update the number of words in that document
        ntokens_in_doc[doc] += count;
        // Add all the tokens
        corpus::token tok; tok.word = word; tok.doc = doc;
        for(size_t i = 0; i < count; ++i) tokens.push_back(tok);
      }
    }
  }
  fin.close();
  ndocs = ntokens_in_doc.size();
  ntokens = tokens.size();
}


void corpus::shuffle_tokens() {
  //  Shuffle the order of the tokens (each word in each document)
  graphlab::random::shuffle(tokens);
}



void split(const corpus& base, corpus& c1, corpus& c2, double prop_c1) {
  std::vector<doc_id_type> doc_perm(base.ndocs);
  for(doc_id_type d = 0; d < doc_perm.size(); ++ d) doc_perm[d] = d;
  graphlab::random::shuffle(doc_perm);
  const size_t maxid(prop_c1 * base.ndocs);
  

  c1.nwords = base.nwords;
  c1.ndocs  = base.ndocs;
  c1.dictionary = base.dictionary;
  c1.ntokens_in_doc = base.ntokens_in_doc;
  
  c2.nwords = base.nwords;
  c2.ndocs  = base.ndocs;
  c2.dictionary = base.dictionary;
  c2.ntokens_in_doc = base.ntokens_in_doc;

  for(size_t i = 0; i < base.tokens.size(); ++i) {
    const corpus::token tok(base.tokens[i]);
    if(doc_perm[tok.doc] < maxid) c1.tokens.push_back(tok);
    else c2.tokens.push_back(tok);
  }
  c1.ntokens = c1.tokens.size();
  c2.ntokens = c2.tokens.size();

}
