#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_SAVE_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_SAVE_HPP
#include <boost/function.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/builtin_parsers.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/fs_util.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/hdfs.hpp>
namespace graphlab {

  namespace graph_ops {
   

  template<typename VertexDataType, typename EdgeDataType, typename Fstream, typename Writer>
  void save_vertex_to_stream(const typename distributed_graph<VertexDataType,EdgeDataType>::vertex_type& vtype,
                              const std::string& targetfilename,
                              Fstream& fout,
                              const Writer &writer) {
    fout << writer.save_vertex(vtype, targetfilename);
  }


  template<typename VertexDataType, typename EdgeDataType, typename Fstream, typename Writer>
  void save_edge_to_stream(const typename distributed_graph<VertexDataType,EdgeDataType>::edge_type& etype,
                              const std::string& targetfilename,
                              Fstream& fout,
                              const Writer &writer) {
    std::string ret = writer.save_edge(etype, targetfilename);
    fout << ret;
  }
  
  template <typename VertexDataType, typename EdgeDataType, typename Writer>
  void save(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
            std::string prefix,
            const Writer &writer,
            bool gzip = true,
            bool save_vertex = true,
            bool save_edge = true,
            size_t files_per_machine = 4) {
    graph.dc().full_barrier();
		if(boost::starts_with(prefix, "hdfs://")) {
      save_to_hdfs(graph, prefix, writer, gzip, save_vertex, save_edge, files_per_machine);
    }
    else {
      save_to_posixfs(graph, prefix, writer, gzip, save_vertex, save_edge, files_per_machine);
    }
    graph.dc().full_barrier();
  }

  template <typename VertexDataType, typename EdgeDataType, typename Writer>
  void save_to_posixfs(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                       std::string prefix,
                       const Writer &writer,
                       bool gzip = true,
                       bool save_vertex = true,
                       bool save_edge = true,
                       size_t files_per_machine = 4) {
    
    typedef distributed_graph<VertexDataType, EdgeDataType> graph_type;
    typedef boost::function<void(typename graph_type::vertex_type)>
                                                        vertex_function_type;
    typedef boost::function<void(typename graph_type::edge_type)>
                                                        edge_function_type;
    typedef std::ofstream base_fstream_type;
    typedef boost::iostreams::filtering_stream<boost::iostreams::output>
                                                        boost_fstream_type;

    // figure out the filenames
    std::vector<std::string> graph_files;
    std::vector<base_fstream_type*> outstreams;
    std::vector<boost_fstream_type*> booststreams;

    graph_files.resize(files_per_machine);
    for(size_t i = 0; i < files_per_machine; ++i) {
      graph_files[i] = prefix + "." + tostr(1 + i + graph.procid() * files_per_machine)
                       + "_of_" + tostr(graph.numprocs() * files_per_machine);
      if (gzip) graph_files[i] += ".gz";
    }

    // create the vector of callbacks

    std::vector<vertex_function_type> vertex_callbacks(graph_files.size());
    std::vector<edge_function_type> edge_callbacks(graph_files.size());
    
    for(size_t i = 0; i < graph_files.size(); ++i) {
      std::cout << "Saving to file: " << graph_files[i] << std::endl;
      // open the stream
      base_fstream_type* out_file = new base_fstream_type(
                                            graph_files[i].c_str(),
                                            std::ios_base::out |
                                                    std::ios_base::binary);
      // attach gzip if the file is gzip
      boost_fstream_type* fout = new boost_fstream_type;
      // Using gzip filter
      if (gzip) fout->push(boost::iostreams::gzip_compressor());
      fout->push(*out_file);

      outstreams.push_back(out_file);
      booststreams.push_back(fout);
      // construct the callback for the parallel for

      vertex_callbacks[i] = boost::bind(
            save_vertex_to_stream<VertexDataType, EdgeDataType,
                                  boost_fstream_type, Writer>,
            _1, graph_files[i], boost::ref(*fout), boost::ref(writer));

      edge_callbacks[i] = boost::bind(
            save_edge_to_stream<VertexDataType, EdgeDataType,
                                boost_fstream_type, Writer>,
            _1, graph_files[i], boost::ref(*fout), boost::ref(writer));
    }

    if (save_vertex) graph.parallel_for_vertices(vertex_callbacks);
    if (save_edge) graph.parallel_for_edges(edge_callbacks);

    // cleanup
    for(size_t i = 0; i < graph_files.size(); ++i) {    
      booststreams[i]->pop();
      if (gzip) booststreams[i]->pop();
      delete booststreams[i];
      delete outstreams[i];
    }
    vertex_callbacks.clear();
    edge_callbacks.clear();
    outstreams.clear();
    booststreams.clear();
  }





  template <typename VertexDataType, typename EdgeDataType, typename Writer>
  void save_to_hdfs(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                    std::string prefix,
                    const Writer &writer,
                    bool gzip = true,
                    bool save_vertex = true,
                    bool save_edge = true,
                    size_t files_per_machine = 4) {

    typedef distributed_graph<VertexDataType, EdgeDataType> graph_type;
    typedef boost::function<void(typename graph_type::vertex_type)>
                                                        vertex_function_type;
    typedef boost::function<void(typename graph_type::edge_type)>
                                                        edge_function_type;
    typedef graphlab::hdfs::fstream base_fstream_type;
    typedef boost::iostreams::filtering_stream<boost::iostreams::output>
                                                        boost_fstream_type;

    // figure out the filenames
    std::vector<std::string> graph_files;
    std::vector<base_fstream_type*> outstreams;
    std::vector<boost_fstream_type*> booststreams;
    
    graph_files.resize(files_per_machine);
    for(size_t i = 0; i < files_per_machine; ++i) {
      graph_files[i] = prefix + "." + tostr(1 + i + graph.procid() * files_per_machine)
                       + "_of_" + tostr(graph.numprocs() * files_per_machine);
      if (gzip) graph_files[i] += ".gz";
    }


    ASSERT_TRUE(hdfs::has_hadoop());
    hdfs& hdfs = hdfs::get_hdfs();
    
    // create the vector of callbacks

    std::vector<vertex_function_type> vertex_callbacks(graph_files.size());
    std::vector<edge_function_type> edge_callbacks(graph_files.size());

    for(size_t i = 0; i < graph_files.size(); ++i) {
      std::cout << "Saving to file: " << graph_files[i] << std::endl;
      // open the stream
      base_fstream_type* out_file = new base_fstream_type(hdfs,
                                                          graph_files[i],
                                                          true);
      // attach gzip if the file is gzip
      boost_fstream_type* fout = new boost_fstream_type;
      // Using gzip filter
      if (gzip) fout->push(boost::iostreams::gzip_compressor());
      fout->push(*out_file);

      outstreams.push_back(out_file);
      booststreams.push_back(fout);
      // construct the callback for the parallel for
      vertex_callbacks[i] = boost::bind(
            save_vertex_to_stream<VertexDataType, EdgeDataType,
                                  boost_fstream_type, Writer>,
            _1, graph_files[i], boost::ref(*fout), boost::ref(writer));

      edge_callbacks[i] = boost::bind(
            save_edge_to_stream<VertexDataType, EdgeDataType,
                                boost_fstream_type, Writer>,
            _1, graph_files[i], boost::ref(*fout), boost::ref(writer));
    }

    if (save_vertex) graph.parallel_for_vertices(vertex_callbacks);
    if (save_edge) graph.parallel_for_edges(edge_callbacks);

    // cleanup
    for(size_t i = 0; i < graph_files.size(); ++i) {
      booststreams[i]->pop();
      if (gzip) booststreams[i]->pop();
      delete booststreams[i];
      delete outstreams[i];
    }
    vertex_callbacks.clear();
    edge_callbacks.clear();
    outstreams.clear();
    booststreams.clear();
  }






  template <typename VertexDataType, typename EdgeDataType>
  void save_structure(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                      std::string prefix, 
                      std::string format,
                      bool gzip = true,
                      size_t files_per_machine = 4) {
    if (format == "snap") {
      save(graph,
           prefix,
           builtin_parsers::tsv_writer<VertexDataType, EdgeDataType>(),
           gzip, false, true, files_per_machine);
    }
    else if (format == "tsv") {
      save(graph,
           prefix,
           builtin_parsers::tsv_writer<VertexDataType, EdgeDataType>(),
           gzip, false, true, files_per_machine);
    }
    else {
      logstream(LOG_ERROR)
            << "Unrecognized Format \"" << format << "\"!" << std::endl;
      return;
    }
  }

} // namespace graph_ops

} // namespace graphlab
#endif
