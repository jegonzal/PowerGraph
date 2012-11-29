#include <graphlab/ui/mongoose/mongoose.h>
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <algorithm>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>
namespace lda {
  typedef int count_type;
  typedef std::vector< count_type > factor_type;
}

namespace std {
  inline lda::factor_type& operator+=(lda::factor_type& lvalue,
      const lda::factor_type& rvalue) {
    if(!rvalue.empty()) {
      if(lvalue.empty()) lvalue = rvalue;
      else {
        for(size_t t = 0; t < lvalue.size(); ++t) lvalue[t] += rvalue[t];
      }
    }
    return lvalue;
  } // end of operator +=
}

namespace lda {

typedef uint16_t topic_id_type;
#define NULL_TOPIC (topic_id_type(-1))
typedef std::vector< topic_id_type > assignment_type;

size_t NTOPICS = 50;
double ALPHA = 0.1;
double BETA = 0.1;
size_t NWORDS = 0;
size_t NDOCS = 0;
size_t NTOKENS = 0;
size_t TOPK = 5;
size_t INTERVAL = 10;
size_t WORDID_OFFSET = 0; 
bool JOIN_ON_ID = true;

float MIMNO_S;
factor_type GLOBAL_TOPIC_COUNT;
std::vector<std::string> DICTIONARY;
boost::unordered_map<std::string, size_t> INVERSE_DICTIONARY;
std::vector<int> LOCKED_WORDS; 

size_t MAX_COUNT = 100;
float BURNIN = -1;


/**
 * \brief Create a token changes event tracker which is reported in
 * the GraphLab metrics dashboard.
 */
DECLARE_EVENT(TOKEN_CHANGES);


struct top_words_type {
  graphlab::mutex lock;
  std::string json_string;
  top_words_type() : 
    json_string("{\n" + json_header_string() + "\tvalues: [] \n }") { }
  inline std::string json_header_string() const {
    return
      "\t\"ntopics\": " + graphlab::tostr(NTOPICS) + ",\n" +
      "\t\"nwords\":  " + graphlab::tostr(NWORDS) + ",\n" +
      "\t\"ndocs\":   " + graphlab::tostr(NDOCS) + ",\n" +
      "\t\"ntokens\": " + graphlab::tostr(NTOKENS) + ",\n" +
      "\t\"alpha\":   " + graphlab::tostr(ALPHA) + ",\n" +
      "\t\"beta\":    " + graphlab::tostr(BETA) + ",\n";
  } // end of json header string
} TOP_WORDS;

std::pair<std::string, std::string>
word_cloud_callback(std::map<std::string, std::string>& varmap) {
  TOP_WORDS.lock.lock();
  const std::pair<std::string, std::string>
    pair("text/html",TOP_WORDS.json_string);
  TOP_WORDS.lock.unlock();
  return pair;
}


void set_alpha(double alphaval) {
  if (alphaval <= 0) ALPHA = 1E-5;
  else ALPHA = alphaval;
}

void set_beta(double betaval) {
  if (betaval <= 0) BETA = 1E-5;
  BETA = betaval;
}


std::pair<std::string, std::string>
set_param_callback(std::map<std::string, std::string>& varmap) {
  std::map<std::string,std::string>::const_iterator iter = varmap.find("alpha");
  if (iter != varmap.end()) {
    double alphaval = atof(iter->second.c_str());
    graphlab::procid_t nprocs = graphlab::dc_impl::get_last_dc()->numprocs();
    graphlab::dc_impl::get_last_dc()->remote_call(boost::counting_iterator<graphlab::procid_t>(0), 
                                        boost::counting_iterator<graphlab::procid_t>(nprocs), 
                                        set_alpha, alphaval);
  }
  iter = varmap.find("beta");
  if (iter != varmap.end()) {
    double betaval = atof(iter->second.c_str());
    graphlab::procid_t nprocs = graphlab::dc_impl::get_last_dc()->numprocs();
    graphlab::dc_impl::get_last_dc()->remote_call(boost::counting_iterator<graphlab::procid_t>(0), 
                                        boost::counting_iterator<graphlab::procid_t>(nprocs), 
                                        set_beta, betaval);
  }
  std::pair<std::string, std::string> pair("text/plain","");
  std::stringstream strm;
  strm << "alpha = " << ALPHA << "\n"
       << "beta = " << BETA << "\n";
  strm.flush();
  pair.second = strm.str();
  return pair;
}


void reset_word_topic_lock() {
  for (size_t i = 0;i < LOCKED_WORDS.size(); ++i) LOCKED_WORDS[i] = -1;
}


void word_topic_lock(size_t wordid, size_t topicid) {
  LOCKED_WORDS[wordid] = topicid;
}

std::pair<std::string, std::string>
lock_word_callback(std::map<std::string, std::string>& varmap) {
  std::pair<std::string, std::string> ret("text/plain","");
  std::map<std::string,std::string>::const_iterator iter = varmap.find("reset");
  if (iter != varmap.end()) {
    ret.second = "reset";
    // reset
    std::cout << "Reset Locks" << std::endl;
    graphlab::procid_t nprocs = graphlab::dc_impl::get_last_dc()->numprocs();
    graphlab::dc_impl::get_last_dc()->remote_call(boost::counting_iterator<graphlab::procid_t>(0), 
                                        boost::counting_iterator<graphlab::procid_t>(nprocs), 
                                        reset_word_topic_lock);
  } else {
    iter = varmap.find("word");
    if (iter != varmap.end()) {
      std::string word = iter->second;
      boost::to_lower(word); 
      boost::trim(word);
      if (INVERSE_DICTIONARY.count(word) == 0) {
        ret.second = "Unable to find word";
        return ret;
      }
      size_t wordid = INVERSE_DICTIONARY[word];
      // get the topic id
      iter = varmap.find("topic");
      size_t topicid = atoi(iter->second.c_str());
      if (topicid >= NTOPICS) {
        ret.second = "Invalid topic number";
        return ret;
      }
      std::cout << "Locking word " << DICTIONARY[wordid] << " to topic " << topicid << std::endl;
      graphlab::procid_t nprocs = graphlab::dc_impl::get_last_dc()->numprocs();
      graphlab::dc_impl::get_last_dc()->remote_call(boost::counting_iterator<graphlab::procid_t>(0), 
                                          boost::counting_iterator<graphlab::procid_t>(nprocs), 
                                          word_topic_lock, wordid, topicid);
      ret.second = "ok";
    }
  }
  return ret;
}


// Graph Types
// ============================================================================

/**
 * \brief The vertex data represents each term and document in the
 * corpus and contains the counts of tokens in each topic.
 */
struct vertex_data {
   size_t join_key; // hash value of the user name, used as vertex join key
  ///! The total number of updates
  uint32_t nupdates;
  ///! The total number of changes to adjacent tokens
  uint32_t nchanges;
  ///! The count of tokens in each topic
  factor_type factor;
  float MIMNO_R;
  vertex_data() : join_key(-1), nupdates(0), nchanges(0), factor(NTOPICS), MIMNO_R(0) { }
  vertex_data(size_t join_key) : join_key(join_key),
                                 nupdates(0), nchanges(0), factor(NTOPICS), MIMNO_R(0) { }
  void save(graphlab::oarchive& arc) const {
    arc << join_key << nupdates << nchanges << MIMNO_R; 
    uint16_t ni = 0;
    for (size_t i = 0;i < factor.size(); ++i) {
      ni += (factor[i] > 0);
    }
    arc << ni;
    for (size_t i = 0;i < factor.size(); ++i) {
      if (factor[i] > 0) {
        arc << uint16_t(i) << factor[i];
      }
    }
  }
  void load(graphlab::iarchive& arc) {
    arc >> join_key >> nupdates >> nchanges >> MIMNO_R; 
    for (size_t i = 0;i < factor.size(); ++i) factor[i] = 0;
    uint16_t ni;
    arc >> ni; 
    for (uint16_t i = 0;i < ni; ++i) {
      uint16_t u; arc >> u;
      arc >> factor[u];
    }
  }
}; // end of vertex_data


/**
 * \brief The edge data represents the individual tokens (word,doc)
 * pairs and their assignment to topics.
 */
struct edge_data {
  ///! The number of changes on the last update
  uint16_t nchanges;
  ///! The assignment of all tokens
  assignment_type assignment;
  edge_data(size_t ntokens = 0) : nchanges(0), assignment(ntokens, NULL_TOPIC) { }
  void save(graphlab::oarchive& arc) const { arc << nchanges << assignment; }
  void load(graphlab::iarchive& arc) { arc >> nchanges >> assignment; }
}; // end of edge_data


/**
 * \brief The LDA graph is a bipartite graph with docs connected to
 * terms if the term occurs in the document.
 *
 * The edges store the number of occurrences of the term in the
 * document as a vector of the assignments of that term in that
 * document to topics.
 *
 * The vertices store the total topic counts.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief The edge_line_parser is used by graph.load to parse lines of the
 * text data file.
 *
 * The global variable MAX_COUNT limits the number of tokens that can
 * be constructed on a particular edge.
 *
 * We use the relativley fast boost::spirit parser to parse each line.
 */
bool edge_line_parser(graph_type& graph, const std::string& fname,
                  const std::string& line) {
  ASSERT_FALSE(line.empty());
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graphlab::vertex_id_type doc_id(-1), word_id(-1);
  size_t count = 0;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(doc_id) = qi::_1] >> -qi::char_(',') >>
      qi::ulong_[phoenix::ref(word_id) = qi::_1] >> -qi::char_(',') >>
      qi::ulong_[phoenix::ref(count) = qi::_1]
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;  
  // Threshold the count
  count = std::min(count, MAX_COUNT);
  // since this is a bipartite graph I need a method to number the
  // left and right vertices differently.  To accomplish I make sure
  // all vertices have non-zero ids and then negate the right vertex.
  // Unfortunatley graphlab reserves -1 and so we add 2 and negate.
  doc_id += 2;
  ASSERT_GT(doc_id, 1);
  doc_id = -doc_id;
  word_id -= WORDID_OFFSET;
  ASSERT_NE(doc_id, word_id);
  // Create an edge and add it to the graph
  graph.add_edge(doc_id, word_id, edge_data(count));
  return true; // successful load
}; // end of graph loader

/*
 * Load a table of userid -> username, used for create vertex_join key
 */
bool vertex_line_parser(graph_type& graph,
                        const std::string& filename,
                        const std::string& line) {
  ASSERT_FALSE(line.empty());
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graphlab::vertex_id_type vid(-1);
  std::string name;

  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(',')  >>  
      +qi::char_[phoenix::ref(name) = qi::_1]
     )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;  

  size_t join_key = boost::hash<std::string>()(name);

  // std::cout << vid << ", " << name << "(" << join_key << ")" << std::endl;

  size_t docid = -(vid+2);
  vertex_data vdata(join_key); 
  ASSERT_FALSE(vdata.join_key == size_t(-1)); // -1 is reserved for non_join_vertex
  graph.add_vertex(docid, vdata);
  return true;
}

inline bool is_word(const graph_type::vertex_type& vertex) {
  return vertex.num_in_edges() > 0 ? 1 : 0;
}

inline bool is_doc(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0 ? 1 : 0;
}

size_t right_emit_key (const graph_type::vertex_type& vertex) {
  return 
    is_doc(vertex) ?
    (JOIN_ON_ID ? (-vertex.id()-2) : vertex.data().join_key) :
    size_t(-1);
}

inline bool is_join(const graph_type::vertex_type& vertex) {
  return right_emit_key(vertex) != size_t(-1);
}

inline size_t count_tokens(const graph_type::edge_type& edge) {
  return edge.data().assignment.size();
}

inline graph_type::vertex_type
get_other_vertex(const graph_type::edge_type& edge,
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}


/**
 <docid> <wordid> <count>
 \verbatim
    0    0     2
    0    4     1
    0    2     3
 \endverbatim
 * 
 * implies that document zero contains word zero twice, word 4 once,
 * and word two three times.
 */
bool load_and_initialize_graph(graphlab::distributed_control& dc,
                               graph_type& graph,
                               const std::string& corpus_dir,
                               const std::string& vid2name_dir
			       ) {
  dc.cout() << "Loading LDA graph." << std::endl;
  graphlab::timer timer; timer.start();

  graph.load(corpus_dir, edge_line_parser);
  graph.load(vid2name_dir, vertex_line_parser);

  dc.cout() << "Loading LDA graph finished in "
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing LDA graph. Finished in "
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Computing number of words and documents." << std::endl;
  NWORDS = graph.map_reduce_vertices<size_t>(is_word);
  NDOCS = graph.map_reduce_vertices<size_t>(is_doc);
  NTOKENS = graph.map_reduce_edges<size_t>(count_tokens);
  size_t isjoin = graph.map_reduce_vertices<size_t>(is_join);
  ASSERT_EQ(isjoin, NDOCS);
  dc.cout() << "Number of words:     " << NWORDS  << std::endl;
  dc.cout() << "Number of docs:      " << NDOCS   << std::endl;
  dc.cout() << "Number of tokens:    " << NTOKENS << std::endl;
  // Prepare the json struct with the word counts
  TOP_WORDS.lock.lock();
  TOP_WORDS.json_string = "{\n" + TOP_WORDS.json_header_string() +
    "\t\"values\": [] \n }";
  TOP_WORDS.lock.unlock();
  return true;
} // end of load and initialize graph


/**
 * \brief Load the dictionary global variable from the file containing
 * the terms (one term per line).
 *
 * Note that while graphs can be loaded from multiple files the
 * dictionary must be in a single file.  The dictionary is loaded
 * entirely into memory and used to display word clouds and the top
 * terms in each topic.
 *
 * \param [in] fname the file containing the dictionary data.  The
 * data can be located on HDFS and can also be gzipped (must end in
 * ".gz").
 * 
 */
bool load_dictionary(const std::string& fname)  {
  // std::cout << "staring load on: "
  //           << graphlab::get_local_ip_as_str() << std::endl;
  const bool gzip = boost::ends_with(fname, ".gz");
  // test to see if the graph_dir is an hadoop path
  if(boost::starts_with(fname, "hdfs://")) {
    graphlab::hdfs hdfs;
    graphlab::hdfs::fstream in_file(hdfs, fname);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    fin.set_auto_close(false);
    if(gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    while(std::getline(fin,term).good()) {
      DICTIONARY.push_back(term);
      INVERSE_DICTIONARY[term] = DICTIONARY.size() - 1;
    }
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } else {
    std::cout << "opening: " << fname << std::endl;
    std::ifstream in_file(fname.c_str(),
                          std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    if (gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good() || !fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    std::cout << "Loooping" << std::endl;
    while(std::getline(fin, term).good()) {
      DICTIONARY.push_back(term);
      INVERSE_DICTIONARY[term] = DICTIONARY.size() - 1;
    }

    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  // std::cout << "Finished load on: "
  //           << graphlab::get_local_ip_as_str() << std::endl;
  std::cout << "Dictionary Size: " << DICTIONARY.size() << std::endl;
  LOCKED_WORDS.resize(DICTIONARY.size(), -1);
  return true;
} // end of load dictionary


// ========================================================
// The Collapsed Gibbs Sampler Function



/**
 * \brief The gather type for the collapsed Gibbs sampler is used to
 * collect the topic counts on adjacent edges so that the apply
 * function can compute the correct topic counts for the center
 * vertex.
 *
 */
struct gather_type {
  factor_type factor;
  uint32_t nchanges;
  gather_type() : nchanges(0) { };
  gather_type(uint32_t nchanges) : factor(NTOPICS), nchanges(nchanges) { };
  void save(graphlab::oarchive& arc) const { arc << factor << nchanges; }
  void load(graphlab::iarchive& arc) { arc >> factor >> nchanges; }
  gather_type& operator+=(const gather_type& other) {
    factor += other.factor;
    nchanges += other.nchanges;
    return *this;
  }
}; // end of gather type







/**
 * \brief The collapsed Gibbs sampler vertex program updates the topic
 * counts for the center vertex and then draws new topic assignments
 * for each edge durring the scatter phase.
 * 
 */
class cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE
   {
public:

  /**
   * \brief At termination we want to disable sampling to allow the
   * correct final counts to be computed.
   */
  static bool DISABLE_SAMPLING; 

   /** \brief gather on all edges */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } // end of gather_edges

  /**
   * \brief Collect the current topic count on each edge.
   */
  gather_type gather(icontext_type& context, const vertex_type& vertex,
                     edge_type& edge) const {
    gather_type ret(edge.data().nchanges);
    const assignment_type& assignment = edge.data().assignment;
    foreach(topic_id_type asg, assignment) {
      if(asg != NULL_TOPIC) ++ret.factor[asg];
    }
    return ret;
  } // end of gather


  /**
   * \brief Update the topic count for the center vertex.  This
   * ensures that the center vertex has the correct topic count before
   * resampling the topics for each token along each edge.
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum) {
    const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
    ASSERT_GT(num_neighbors, 0);
    // There should be no new edge data since the vertex program has been cleared
    vertex_data& vdata = vertex.data();
    ASSERT_EQ(sum.factor.size(), NTOPICS);
    ASSERT_EQ(vdata.factor.size(), NTOPICS);
    vdata.nupdates++;
    vdata.nchanges = sum.nchanges;
    vdata.factor = sum.factor;
    if (is_doc(vertex)) {
      float MIMNO_R = 0.0;
      for (size_t i = 0;i < vdata.factor.size(); ++i) {
        MIMNO_R += vdata.factor[i] * BETA / (BETA * NWORDS + GLOBAL_TOPIC_COUNT[i]);
      }
      vdata.MIMNO_R = MIMNO_R;
    }
 } // end of apply


  /**
   * \brief Scatter on all edges if the computation is on-going.
   * Computation stops after bunrin or when disable sampling is set to
   * true.
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return (DISABLE_SAMPLING || (BURNIN > 0 && context.elapsed_seconds() > BURNIN))? 
      graphlab::NO_EDGES : graphlab::ALL_EDGES;
  }; // end of scatter edges


  /**
   * \brief Draw new topic assignments for each edge token.
   *
   * Note that we exploit the GraphLab caching model here by DIRECTLY
   * modifying the topic counts of adjacent vertices.  Making the
   * changes immediately visible to any adjacent vertex programs
   * running on the same machine.  However, these changes will be
   * overwritten during the apply step and are only used to accelerate
   * sampling.  This is a potentially dangerous violation of the
   * abstraction and should be taken with caution.  In our case all
   * vertex topic counts are preallocated and atomic operations are
   * used.  In addition during the sampling phase we must be careful
   * to guard against potentially negative temporary counts.
   */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    factor_type& doc_topic_count =  is_doc(edge.source()) ?
      edge.source().data().factor : edge.target().data().factor;
    factor_type& word_topic_count = is_word(edge.source()) ?
      edge.source().data().factor : edge.target().data().factor;
    ASSERT_EQ(doc_topic_count.size(), NTOPICS);
    ASSERT_EQ(word_topic_count.size(), NTOPICS);
    float MIMNO_R = is_doc(edge.source()) ? edge.source().data().MIMNO_R :
                      edge.target().data().MIMNO_R;
    float MIMNO_Q = 0.0;
    std::vector<float> MIMNO_Q_CACHE(NTOPICS);
    
    size_t wordid = is_word(edge.source()) ? edge.source().id() : edge.target().id();

    for (size_t t = 0; t < NTOPICS; ++t) {
      const float n_wt  =
        std::max(count_type(word_topic_count[t]), count_type(0));
     if (n_wt > 0) {
      const float n_dt =
          std::max(count_type(doc_topic_count[t]), count_type(0));
      const float n_t  =
        std::max(count_type(GLOBAL_TOPIC_COUNT[t]), count_type(0));
       MIMNO_Q_CACHE[t] = (ALPHA + n_dt)/(BETA * NWORDS + n_t); 
       MIMNO_Q_CACHE[t] = MIMNO_Q_CACHE[t] * n_wt; 
       MIMNO_Q += MIMNO_Q_CACHE[t]; 
     }
    }

    // run the actual gibbs sampling
    std::vector<float> prob(NTOPICS);
    assignment_type& assignment = edge.data().assignment;
    edge.data().nchanges = 0;

    foreach(topic_id_type& asg, assignment) {
      const topic_id_type old_asg = asg;
      if(asg != NULL_TOPIC) { // construct the cavity
        --doc_topic_count[asg];
        --word_topic_count[asg];
        --GLOBAL_TOPIC_COUNT[asg];
      const float n_dt =
          std::max(count_type(doc_topic_count[asg]), count_type(0));
      const float n_t  =
        std::max(count_type(GLOBAL_TOPIC_COUNT[asg]), count_type(0));
      const float n_wt  =
        std::max(count_type(word_topic_count[asg]), count_type(0));


        MIMNO_Q -= MIMNO_Q_CACHE[asg];
        MIMNO_Q_CACHE[asg] = (ALPHA + n_dt)/(BETA * NWORDS + n_t) * n_wt;
        MIMNO_Q += MIMNO_Q_CACHE[asg]; 
      }
      asg = 0; 
      ASSERT_GE(MIMNO_S, 0);
      ASSERT_GE(MIMNO_R, 0);
      ASSERT_GE(MIMNO_Q, 0);
      float f = graphlab::random::uniform<float>(0, MIMNO_S + MIMNO_R + MIMNO_Q);
      if (f < MIMNO_S) {
        float ctr = 0;
        
        for (size_t t = 0; t < NTOPICS; ++t) {
          ctr += ALPHA * BETA / (BETA * NWORDS + GLOBAL_TOPIC_COUNT[t]);
          if (ctr >= f) {
            asg = t;
            break;
          }
        }
      }
      else if (f < MIMNO_S + MIMNO_R) {
        float ctr = 0;
        f = f - MIMNO_S;
        for(size_t t = 0; t < NTOPICS; ++t) {
          if (doc_topic_count[t] > 0) {
            ctr += doc_topic_count[t] * BETA / (BETA * NWORDS + GLOBAL_TOPIC_COUNT[t]);
            if (ctr >= f) {
              asg = t;
              break;
            }
          }
        }
      }
      else {
        f = f - MIMNO_S - MIMNO_R;
        float ctr = 0;
        for(size_t t = 0; t < NTOPICS; ++t) {
          if (word_topic_count[t] > 0) {
            ctr += MIMNO_Q_CACHE[t]; 
            if (ctr >= f) {
              asg = t;
              break;
            }
          }
        }
      }
      if (LOCKED_WORDS[wordid] != -1) asg = LOCKED_WORDS[wordid];
      // asg = std::max_element(prob.begin(), prob.end()) - prob.begin();
      ++doc_topic_count[asg];
      ++word_topic_count[asg];
      ++GLOBAL_TOPIC_COUNT[asg];
      MIMNO_Q -= MIMNO_Q_CACHE[asg];
{
      const float n_dt =
          std::max(count_type(doc_topic_count[asg]), count_type(0));
      const float n_t  =
        std::max(count_type(GLOBAL_TOPIC_COUNT[asg]), count_type(0));
      const float n_wt  =
        std::max(count_type(word_topic_count[asg]), count_type(0));

      MIMNO_Q_CACHE[asg] = (ALPHA + n_dt)/(BETA * NWORDS + n_t) * n_wt;
      MIMNO_Q += MIMNO_Q_CACHE[asg]; 
}
      if(asg != old_asg) {
        ++edge.data().nchanges;
      }
      INCREMENT_EVENT(TOKEN_CHANGES,1);
    } // End of loop over each token
    // singla the other vertex
    context.signal(get_other_vertex(edge, vertex));
  } // end of scatter function

}; // end of cgs_lda_vertex_program


bool cgs_lda_vertex_program::DISABLE_SAMPLING = false;


/**
 * \brief The icontext type associated with the cgs_lda_vertex program
 * is needed for all aggregators.
 */
typedef cgs_lda_vertex_program::icontext_type icontext_type;


// ========================================================
// Aggregators


/**
 * \brief The topk aggregator is used to periodically compute and
 * display the topk most common words in each topic.
 *
 * The number of words is determined by the global variable \ref TOPK
 * and the interval is determined by the global variable \ref INTERVAL.
 *
 */
class topk_aggregator {
  typedef std::pair<float, graphlab::vertex_id_type> cw_pair_type;
private:
  std::vector< std::set<cw_pair_type> > top_words;
  size_t nchanges, nupdates;
public:
  topk_aggregator(size_t nchanges = 0, size_t nupdates = 0) :
    nchanges(nchanges), nupdates(nupdates) { }

  void save(graphlab::oarchive& arc) const { arc << top_words << nchanges; }
  void load(graphlab::iarchive& arc) { arc >> top_words >> nchanges; }

  topk_aggregator& operator+=(const topk_aggregator& other) {
    nchanges += other.nchanges;
    nupdates += other.nupdates;
    if(other.top_words.empty()) return *this;
    if(top_words.empty()) top_words.resize(NTOPICS);
    for(size_t i = 0; i < top_words.size(); ++i) {
      // Merge the topk
      top_words[i].insert(other.top_words[i].begin(),
                          other.top_words[i].end());
      // Remove excess elements
      while(top_words[i].size() > TOPK)
        top_words[i].erase(top_words[i].begin());
    }
    return *this;
  } // end of operator +=

  static topk_aggregator map(icontext_type& context,
                             const graph_type::vertex_type& vertex) {
    topk_aggregator ret_value;
    const vertex_data& vdata = vertex.data();
    ret_value.nchanges = vdata.nchanges;
    ret_value.nupdates = vdata.nupdates;
    if(is_word(vertex)) {
      const graphlab::vertex_id_type wordid = vertex.id();
      ret_value.top_words.resize(vdata.factor.size());
      for(size_t i = 0; i < vdata.factor.size(); ++i) {
        const cw_pair_type pair(vdata.factor[i], wordid);
        ret_value.top_words[i].insert(pair);
      }
    }
    return ret_value;
  } // end of map function


  static void finalize(icontext_type& context,
                       const topk_aggregator& total) {
    if(context.procid() != 0) return;
    std::string json = "{\n"+ TOP_WORDS.json_header_string() +
      "\t\"values\": [\n";
    for(size_t i = 0; i < total.top_words.size(); ++i) {
      std::cout << "Topic " << i << ": ";
      json += "\t[\n";
      size_t counter = 0;
      rev_foreach(cw_pair_type pair, total.top_words[i])  {
      ASSERT_LT(pair.second, DICTIONARY.size());
        json += "\t\t[\"" + DICTIONARY[pair.second] + "\", " +
          graphlab::tostr(pair.first) + "]";
        if(++counter < total.top_words[i].size()) json += ", ";
        json += '\n';
        std::cout << DICTIONARY[pair.second]
                  << "(" << pair.first << ")" << ", ";
        // std::cout << DICTIONARY[pair.second] << ",  ";
      }
      json += "\t]";
      if(i+1 < total.top_words.size()) json += ", ";
      json += '\n';
      std::cout << std::endl;
    }
    json += "]}";
    // Post the change to the global variable
    TOP_WORDS.lock.lock();
    TOP_WORDS.json_string.swap(json);
    TOP_WORDS.lock.unlock();

    std::cout << "\nNumber of token changes: " << total.nchanges << std::endl;
    std::cout << "\nNumber of updates:       " << total.nupdates << std::endl;
  } // end of finalize
}; // end of topk_aggregator struct



/**
 * \brief The global counts aggregator computes the total number of
 * tokens in each topic across all words and documents and then
 * updates the \ref GLOBAL_TOPIC_COUNT variable.
 *
 */
struct global_counts_aggregator {
  typedef graph_type::vertex_type vertex_type;
  static factor_type map(icontext_type& context, const vertex_type& vertex) {
    return vertex.data().factor;
  } // end of map function

  static void finalize(icontext_type& context, const factor_type& total) {
    size_t sum = 0;
    float NEW_MIMNO_S = 0;
    for(size_t t = 0; t < total.size(); ++t) {
      GLOBAL_TOPIC_COUNT[t] =
        std::max(count_type(total[t]/2), count_type(0));
      sum += GLOBAL_TOPIC_COUNT[t];
      NEW_MIMNO_S += ALPHA * BETA / (BETA * NWORDS + (GLOBAL_TOPIC_COUNT[t] > 0 ? GLOBAL_TOPIC_COUNT[t] : 0));
    }
    MIMNO_S = NEW_MIMNO_S;
    context.cout() << "Total Tokens: " << sum << std::endl;
  } // end of finalize
}; // end of global_counts_aggregator struct




/**
 * \brief The Likelihood aggregators maintains the current estimate of
 * the log-likelihood of the current token assignments.
 *
 *  llik_words_given_topics = ...
 *    ntopics * (gammaln(nwords * beta) - nwords * gammaln(beta)) - ...
 *    sum_t(gammaln( n_t + nwords * beta)) +
 *    sum_w(sum_t(gammaln(n_wt + beta)));
 *
 *  llik_topics = ...
 *    ndocs * (gammaln(ntopics * alpha) - ntopics * gammaln(alpha)) + ...
 *    sum_d(sum_t(gammaln(n_td + alpha)) - gammaln(sum_t(n_td) + ntopics * alpha));
 */
class likelihood_aggregator : public graphlab::IS_POD_TYPE {
  typedef graph_type::vertex_type vertex_type;
  double lik_words_given_topics;
  double lik_topics;
public:
  likelihood_aggregator() : lik_words_given_topics(0), lik_topics(0) { }

  likelihood_aggregator& operator+=(const likelihood_aggregator& other) {
    lik_words_given_topics += other.lik_words_given_topics;
    lik_topics += other.lik_topics;
    return *this;
  } // end of operator +=

  static likelihood_aggregator
  map(icontext_type& context, const vertex_type& vertex) {
    using boost::math::lgamma;
    const factor_type& factor = vertex.data().factor;
    ASSERT_EQ(factor.size(), NTOPICS);
   likelihood_aggregator ret;
    if(is_word(vertex)) {
      for(size_t t = 0; t < NTOPICS; ++t) {
        const double value = std::max(count_type(factor[t]), count_type(0));
        ret.lik_words_given_topics += lgamma(value + BETA);
      }
    } else {  ASSERT_TRUE(is_doc(vertex));
      double ntokens_in_doc = 0;
      for(size_t t = 0; t < NTOPICS; ++t) {
        const double value = std::max(count_type(factor[t]), count_type(0));
        ret.lik_topics += lgamma(value + ALPHA);
        ntokens_in_doc += factor[t];
      }
      ret.lik_topics -= lgamma(ntokens_in_doc + NTOPICS * ALPHA);
    }
    return ret;
  } // end of map function

  static void finalize(icontext_type& context, const likelihood_aggregator& total) {
    using boost::math::lgamma;
    // Address the global sum terms
    double denominator = 0;
    for(size_t t = 0; t < NTOPICS; ++t) {
      denominator += lgamma(GLOBAL_TOPIC_COUNT[t] + NWORDS * BETA);
    } // end of for loop

    const double lik_words_given_topics =
      NTOPICS * (lgamma(NWORDS * BETA) - NWORDS * lgamma(BETA)) -
      denominator + total.lik_words_given_topics;

    const double lik_topics =
      NDOCS * (lgamma(NTOPICS * ALPHA) - NTOPICS * lgamma(ALPHA)) +
      total.lik_topics;

    const double lik = lik_words_given_topics + lik_topics;
    context.cout() << "Likelihood: " << lik << std::endl;
  } // end of finalize
}; // end of likelihood_aggregator struct



/**
 * \brief The selective signal functions are used to signal only the
 * vertices corresponding to words or documents.  This is done by
 * using the iengine::map_reduce_vertices function.
 */
struct signal_only {
  /**
   * \brief Signal only the document vertices and skip the word
   * vertices.
   */ 
  static graphlab::empty
  docs(icontext_type& context, const graph_type::vertex_type& vertex) {
    if(is_doc(vertex)) context.signal(vertex);
    return graphlab::empty();
  } // end of signal_docs
 
 /**
  * \brief Signal only the word vertices and skip the document
  * vertices.
  */
  static graphlab::empty
  words(icontext_type& context, const graph_type::vertex_type& vertex) {
    if(is_word(vertex)) context.signal(vertex);
    return graphlab::empty();
  } // end of signal_words
}; // end of selective_only



struct count_saver {
  bool save_words;
  count_saver(bool save_words) : save_words(save_words) { }
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    // Skip saving vertex data if the vertex type is not consistent
    // with the save type
    if((save_words && is_doc(vertex)) ||
       (!save_words && is_word(vertex))) return "";
    // Proceed to save
    std::stringstream strm;
    if(save_words) {
      const graphlab::vertex_id_type vid = vertex.id();
      strm << vid << '\t';
    } else { // save documents
      const graphlab::vertex_id_type vid = (-vertex.id()) - 2;
      strm << vid << '\t';
    }
    const factor_type& factor = vertex.data().factor;
    for(size_t i = 0; i < factor.size(); ++i) { 
      strm << factor[i];
      if(i+1 < factor.size()) strm << '\t';
    }
    strm << '\n';
    return strm.str();
  }
  std::string save_edge(const edge_type& edge) const {
    return ""; //nop
  }
}; // end of prediction_saver


/**
 * \brief The omni engine type is used to allow switching between
 * synchronous and asynchronous computation. 
 */
typedef graphlab::omni_engine<cgs_lda_vertex_program> engine_type;

void initialize_global() {
  ADD_CUMULATIVE_EVENT(TOKEN_CHANGES, "Token Changes", "Changes");
  lda::GLOBAL_TOPIC_COUNT.resize(NTOPICS);
}
// 
// 
// int main(int argc, char** argv) {
//   global_logger().set_log_level(LOG_INFO);
//   global_logger().set_log_to_console(true);
//   ///! Initialize control plain using mpi
//   graphlab::mpi_tools::init(argc, argv);
//   graphlab::distributed_control dc;
//   //  INITIALIZE_EVENT_LOG(dc);
//   ADD_CUMULATIVE_EVENT(TOKEN_CHANGES, "Token Changes", "Changes");
// 
//   graphlab::command_line_options clopts(description);
//   std::string corpus_dir;
//   std::string dictionary_fname;
//   std::string doc_dir;
//   std::string word_dir;
//   std::string exec_type = "asynchronous";
//   std::string format = "matrix";
//   clopts.attach_option("dictionary", dictionary_fname,
//                        "The file containing the list of unique words");
//   clopts.attach_option("engine", exec_type, 
//                        "The engine type synchronous or asynchronous");
//   clopts.attach_option("corpus", corpus_dir,
//                        "The directory or file containing the corpus data.");
//   clopts.add_positional("corpus");
//   clopts.attach_option("ntopics", NTOPICS,
//                        "Number of topics to use.");
//   clopts.attach_option("alpha", ALPHA,
//                        "The document hyper-prior");
//   clopts.attach_option("beta", BETA,
//                        "The word hyper-prior");
//   clopts.attach_option("topk", TOPK,
//                        "The number of words to report");
//   clopts.attach_option("interval", INTERVAL,
//                        "statistics reporting interval");
//   clopts.attach_option("max_count", MAX_COUNT,
//                        "The maximum number of occurences of a word in a document.");
//   clopts.attach_option("format", format,
//                        "Formats: matrix,json,json-gzip");
//   clopts.attach_option("burnin", BURNIN, 
//                        "The time in second to run until a sample is collected. "
//                        "If less than zero the sampler runs indefinitely.");
//   clopts.attach_option("doc_dir", doc_dir,
//                        "The output directory to save the final document counts.");
//   clopts.attach_option("word_dir", word_dir,
//                        "The output directory to save the final words counts.");
// 
// 
//   if(!clopts.parse(argc, argv)) {
//     graphlab::mpi_tools::finalize();
//     return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
//   }
// 
//   if(dictionary_fname.empty()) {
//     logstream(LOG_WARNING) << "No dictionary file was provided." << std::endl
//                            << "Top k words will not be estimated." << std::endl;
//   }
// 
//   if(corpus_dir.empty()) {
//     logstream(LOG_ERROR) << "No corpus file was provided." << std::endl;
//     return EXIT_FAILURE;
//   }
// 
//   
//   ///! Initialize global variables
//   GLOBAL_TOPIC_COUNT.resize(NTOPICS);
//   if(!dictionary_fname.empty()) {
//     const bool success = load_dictionary(dictionary_fname);
//     if(!success) {
//       logstream(LOG_ERROR) << "Error loading dictionary." << std::endl;
//       return EXIT_FAILURE;
//     }
//   }
// 
//   ///! load the graph
//   graph_type graph(dc, clopts);
//   {
//     const bool success = 
//       load_and_initialize_graph(dc, graph, corpus_dir, format);
//     if(!success) {
//       logstream(LOG_ERROR) << "Error loading graph." << std::endl;
//       return EXIT_FAILURE;
//     }
//   }
// 
// 
//   const size_t ntokens = graph.map_reduce_edges<size_t>(count_tokens);
//   dc.cout() << "Total tokens: " << ntokens << std::endl;
// 
// 
// 
//   engine_type engine(dc, graph, exec_type, clopts);
//   ///! Add an aggregator
//   if(!DICTIONARY.empty()) {
//     const bool success =
//       engine.add_vertex_aggregator<topk_aggregator>
//       ("topk", topk_aggregator::map, topk_aggregator::finalize) &&
//       engine.aggregate_periodic("topk", INTERVAL);
//     ASSERT_TRUE(success);
//   }
// 
//   { // Add the Global counts aggregator
//     const bool success =
//       engine.add_vertex_aggregator<factor_type>
//       ("global_counts", 
//        global_counts_aggregator::map, 
//        global_counts_aggregator::finalize) &&
//       engine.aggregate_periodic("global_counts", 5);
//     ASSERT_TRUE(success);
//   }
//   
// /*  { // Add the likelihood aggregator
//     const bool success =
//       engine.add_vertex_aggregator<likelihood_aggregator>
//       ("likelihood", 
//        likelihood_aggregator::map, 
//        likelihood_aggregator::finalize) &&
//       engine.aggregate_periodic("likelihood", 10);
//     ASSERT_TRUE(success);
//   }*/
// 
//   ///! schedule only documents
//   dc.cout() << "Running The Collapsed Gibbs Sampler" << std::endl;
//   engine.map_reduce_vertices<graphlab::empty>(signal_only::docs);
//   graphlab::timer timer;
//   // Enable sampling
//   cgs_lda_vertex_program::DISABLE_SAMPLING = false;
//   // Run the engine
//   engine.start();
//   // Finalize the counts
//   cgs_lda_vertex_program::DISABLE_SAMPLING = true;
//   engine.signal_all();
//   engine.start();
//   
//   const double runtime = timer.current_time();
//   dc.cout()
//     << "----------------------------------------------------------" << std::endl
//     << "Final Runtime (seconds):   " << runtime
//     << std::endl
//     << "Updates executed: " << engine.num_updates() << std::endl
//     << "Update Rate (updates/second): "
//     << engine.num_updates() / runtime << std::endl;
//   
//   
//   
//   if(!word_dir.empty()) {
//     // save word topic counts
//     const bool gzip_output = false;
//     const bool save_vertices = true;
//     const bool save_edges = false;
//     const size_t threads_per_machine = 2;
//     const bool save_words = true;
//     graph.save(word_dir, count_saver(save_words),
//                gzip_output, save_vertices, 
//                save_edges, threads_per_machine);
//   }
// 
//   
//   if(!doc_dir.empty()) {
//     // save doc topic counts
//     const bool gzip_output = false;
//     const bool save_vertices = true;
//     const bool save_edges = false;
//     const size_t threads_per_machine = 2;
//     const bool save_words = false;
//     graph.save(doc_dir, count_saver(save_words),
//                gzip_output, save_vertices, 
//                save_edges, threads_per_machine);
// 
//   }
// 
// 
//   graphlab::stop_metric_server_on_eof();
//   graphlab::mpi_tools::finalize();
//   return EXIT_SUCCESS;
// } // end of main
} // end of namespace
#include <graphlab/macros_undef.hpp>
