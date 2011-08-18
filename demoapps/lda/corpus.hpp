#ifndef CORPUS_HPP
#define CORPUS_HPP

#include <stdint.h>
#include <cassert>
#include <string>
#include <vector>


typedef uint32_t word_id_type;
typedef uint32_t doc_id_type;
typedef uint16_t topic_id_type;

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
