#ifndef CORPUS
#define CORPUS

#include <stdint.h>
#include <cassert>
#include <string>
#include <vector>


typedef uint32_t word_id_type;
typedef uint32_t doc_id_type;
typedef uint16_t topic_id_type;

struct token_type {
  word_id_type word;
  doc_id_type doc;
  token_type(const word_id_type& word = 0, const doc_id_type& doc = 0) : 
    word(word), doc(doc) { }
};

std::ostream& operator<<(std::ostream& out, const token_type& tok);

#define NULL_TOPIC topic_id_type(-1)


struct corpus_type {
  size_t nwords;
  size_t ndocs;
  size_t ntokens;


  std::vector< token_type > tokens;
  std::vector<std::string> dictionary;
  std::vector< word_id_type > ntokens_in_doc;

  corpus_type(const std::string& dictionary_fname, 
              const std::string& counts_fname );

  void load_dictionary(const std::string& fname);

  void load_counts(const std::string& fname);

  void shuffle_tokens();

};

/**
 * Randomly split the corpus into two separate sets
 */
void split(const corpus_type& base, corpus_type& c1, corpus_type& c2, double prop_c1);



#endif
