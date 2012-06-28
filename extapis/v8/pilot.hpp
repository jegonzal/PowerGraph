#ifndef GRAPHLAB_PILOT_HPP
#define GRAPHLAB_PILOT_HPP

#include <v8.h>
#include <graphlab.hpp>

namespace graphlab {

  class templates;

  // TODO: refactor graph operations to another graph object
  /**
   * Script driver (exposed javascript interface.)
   *
   * Example:
   *    var x = new pilot();
   *    x.ping();     //=> pong
   */
  class pilot {
  public:
    
    // TODO: allow more generic types
    typedef double vertex_data_type;
    typedef double edge_data_type;
    typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
    typedef double gather_type;

  private:
   
    static graphlab_options opts;   // TODO: this is probably wrong
    static templates templs; 

    distributed_control dc;
    graph_type graph;

  public:
 
    pilot();

    /** Prints "pong" to STDOUT */
    void ping();

    /**
     * Loads a graph from the specified path.
     * @param path      path to graph file
     * @param format    format of graph file
     */
    void load_graph(const std::string &path, const std::string &format);

    /**
     * Takes a javascript constructor to create vertex programs.
     */
    void fly(const v8::Handle<v8::Function> &function);

    /**
     * Maps the given javascript function across all vertices.
     */
    void transform_vertices(const v8::Handle<v8::Function> &function);

    /**
     * Adds a JS binding of the class to the given object. Throws
     * a native exception on error.
     * @param dest      v8 global object to bind to
     */
    static void setup_bindings(const v8::Handle<v8::Object> &dest);

    /**
     * Saves command line options for this session.
     * @param clopts    command line options
     * TODO there has got to be a better way than this
     */
    static void set_clopts(const graphlab_options &clopts);

    static templates &get_templates();

  };

  class js_proxy :
  public ivertex_program<pilot::graph_type, pilot::gather_type>,
  public IS_POD_TYPE {
    
    double last_change;
    v8::Persistent<v8::Object> jsobj;
  
  public:  
    
    // TODO: this is a hack
    static v8::Persistent<v8::Function> constructor;
    static void set_ctor(const v8::Handle<v8::Function> &ctor);

    js_proxy();
    pilot::gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const; 
    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total);
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const;
    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const;
    // TODO: other vertex program stuff
    // TODO: pass the context along!

  };

  /**
   * Wrapper for a JS function that is suitable for passing to
   * transform_vertices, transform_edges etc.
   */
  struct js_functor {
    static v8::Persistent<v8::Function> function;
    // not really a functor
    static void set_function(const v8::Handle<v8::Function>& func);
    static void invoke(pilot::graph_type::vertex_type& vertex);
  };

};

namespace cvv8 {

  ////////////////////////// V8 POLICIES ///////////////////////////

  // A helper to support converting from pilot to its JS handle.
  typedef NativeToJSMap<graphlab::pilot> BMap;

  // pilot Ctors we want to bind to v8 
  typedef Signature<graphlab::pilot (CtorForwarder<graphlab::pilot *()>)> pilotCtors;

  CVV8_TypeName_DECL((graphlab::pilot));

 /**
  * This policy class is required unless you just want to bind to 
  * the default constructor. It creates native objects for the 
  * underlying binding code.
  * @internal
  */
  template <>
  class ClassCreator_Factory<graphlab::pilot> {
  public:
      typedef graphlab::pilot * ReturnType;
      static ReturnType Create( v8::Persistent<v8::Object> & jsSelf, v8::Arguments const & argv );
      static void Delete( ReturnType obj );
  };

 /**
  * Disable subclassing from within javascript
  * @internal
  */
  template <>
  struct ClassCreator_SearchPrototypeForThis<graphlab::pilot> : Opt_Bool<false> {};
    
 /**
  * Required specialization so that the conversion API can derive
  * the native 'this' object from v8::Arguments::This() and from
  * function arguments of our bound type.
  *
  * This implementation works by using the plumbing installed
  * by ClassCreator.
  * @internal
  */
  template <>
  struct JSToNative<graphlab::pilot>
    : JSToNative_ClassCreator<graphlab::pilot> {};

 /**
  * Native-to-JS conversion. This conversion is only possible when
  * we explicitly add support to the class-binding code to add
  * the necessary binding metadata. Alternately, if the bound class
  * contains v8-defined data types, e.g. a Handle<Object> referring
  * to itself then implementing NativeToJS is easy to do - just
  * return the handle held by the native.
  *   
  * In this example we're using the NativeToJSMap helper code to 
  * plug in/unplug our bindings during native object 
  * construction/destruction via the the ClassCreator's factory 
  * policy (ClassCreator_Factory<graphlab::pilot>).        
  * @internal
  */
  template <>
  struct NativeToJS<graphlab::pilot>
      :  NativeToJSMap<graphlab::pilot>::NativeToJSImpl {};

};

#endif
