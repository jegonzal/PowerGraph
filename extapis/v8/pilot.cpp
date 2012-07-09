#include "templates.hpp"
#include "pilot.hpp"

using namespace graphlab;
using namespace v8;
namespace cv = cvv8;

/////////////////////////////////////////////////////////////////////////////
// pilot
/////////////////////////////////////////////////////////////////////////////

pilot::pilot() : graph(*dc, opts) {}

void pilot::ping(){ std::cout << "pong." << std::endl; }

void pilot::load_graph(const std::string &path, const std::string &format){
  graph.load_format(path, format);
  graph.finalize();
}

void pilot::load_synthetic_powerlaw(size_t powerlaw){
  graph.load_synthetic_powerlaw(powerlaw);
  graph.finalize();
}

void pilot::fly(const Handle<Function> &function){

  js_proxy::set_ctor(function); // FIXME: should not be a static!
  omni_engine<js_proxy> engine(*dc, graph, "synchronous", opts);
  engine.signal_all(); // TODO: allow user to specify an array of vertices to signal, or all
  engine.start();

  logstream(LOG_INFO) << "done." << std::endl;

  // TODO: move to another function
  const float runtime = engine.elapsed_seconds();
  size_t update_count = engine.num_updates();
  logstream(LOG_EMPH) << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second." << std::endl
            << "Iterations: " << engine.iteration() << std::endl;

}

void pilot::transform_vertices(const Handle<Function> &function){
  js_functor::set_function(function); // FIXME: should not be static!
  graph.transform_vertices(js_functor::invoke);
}

void pilot::save_graph(const std::string &prefix,
                       const Handle<Function> &vwriter,
                       const Handle<Function> &ewriter, // TODO: get rid of booleans
                       bool gzip, bool save_vertex, bool save_edge){
  js_writer writer(vwriter, ewriter);
  graph.save(prefix, writer, gzip, save_vertex, save_edge);
}

/////////////////////////////////////////////////////////////////////////////
// pilot::static
/////////////////////////////////////////////////////////////////////////////

distributed_control *pilot::dc;
graphlab_options pilot::opts;
templates pilot::templs;

void pilot::setup_bindings(const Handle<Object> &dest){
  cv::ClassCreator<pilot>::Instance().SetupBindings(dest);
}

void pilot::set_clopts(const graphlab_options &clopts){ opts = clopts; }
void pilot::set_dc(distributed_control &control){ dc = &control; }

templates &pilot::get_templates(){ return templs; }

/////////////////////////////////////////////////////////////////////////////
// js_proxy
/////////////////////////////////////////////////////////////////////////////

js_proxy::js_proxy() {
  // TODO deal with multi-threaded environments
  HandleScope handle_scope;
  jsobj = Persistent<Object>::New(constructor->NewInstance());
}

js_proxy::js_proxy(const js_proxy& other){
  HandleScope handle_scope;
  this->jsobj = Persistent<Object>::New(other.jsobj);
}

js_proxy &js_proxy::operator=(const js_proxy& other){
  HandleScope handle_scope;
  if (this == &other) return *this;
  this->jsobj.Dispose();
  this->jsobj = Persistent<Object>::New(other.jsobj);
  return *this;
}

js_proxy::~js_proxy(){
  jsobj.Dispose();
}

pilot::gather_type js_proxy::gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
  HandleScope handle_scope;
  Local<Function> f = Function::Cast(*jsobj->Get(JSTR("gather")));
  Handle<Value> ret = cv::CallForwarder<2>::Call(jsobj, f, vertex, edge);
  return cv::CastFromJS<gather_type>(ret);
}

void js_proxy::apply(icontext_type& context, vertex_type& vertex, const gather_type& total){
  HandleScope handle_scope;
  Local<Function> f = Function::Cast(*jsobj->Get(JSTR("apply")));
  cv::CallForwarder<2>::Call(jsobj, f, vertex, total);
}

edge_dir_type js_proxy::scatter_edges(icontext_type& context, const vertex_type& vertex) const {
  HandleScope handle_scope;
  Local<Function> f = Function::Cast(*jsobj->Get(JSTR("scatter_edges")));
  return cv::JSToNative<edge_dir_type>()(cv::CallForwarder<1>::Call(jsobj, f, vertex));
}

void js_proxy::scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
  HandleScope handle_scope;
  Local<Function> f = Function::Cast(*jsobj->Get(JSTR("scatter")));
  Handle<Value> args[] = {
    cv::CastToJS(context), cv::CastToJS(vertex), cv::CastToJS(edge)
  };
  f->Call(jsobj, 3, args);
}

/////////////////////////////////////////////////////////////////////////////
// js_proxy::static
/////////////////////////////////////////////////////////////////////////////

Persistent<Function> js_proxy::constructor;

void js_proxy::set_ctor(const Handle<Function> &ctor){
  // TODO: worry about memory management (should dispose?)
  HandleScope handle_scope;
  constructor  = Persistent<Function>::New(ctor);
}

/////////////////////////////////////////////////////////////////////////////
// js_functor
/////////////////////////////////////////////////////////////////////////////

// TODO: fix this
Persistent<Function> js_functor::function;

void js_functor::invoke(pilot::graph_type::vertex_type &vertex){
  cv::CallForwarder<1>::Call(function, vertex);
}

void js_functor::set_function(const Handle<Function> &func){
  // TODO: worry about memory management
  function = Persistent<Function>::New(func);
}

js_writer::js_writer(const Handle<Function> &vwriter, const Handle<Function> &ewriter){
  vertex_writer = Persistent<Function>::New(vwriter);
  edge_writer = Persistent<Function>::New(ewriter);
}

js_writer::~js_writer(){
  vertex_writer.Dispose();
  edge_writer.Dispose();
}

std::string js_writer::save_vertex(cv::vertex_type vertex){
  return cv::CastFromJS<std::string>(cv::CallForwarder<1>::Call(vertex_writer, vertex));
}

std::string js_writer::save_edge(cv::edge_type edge){
  return cv::CastFromJS<std::string>(cv::CallForwarder<1>::Call(edge_writer, edge));
}

/////////////////////////////////////////////////////////////////////////////
// JS Bindings and Policies
/////////////////////////////////////////////////////////////////////////////

static void bind_member_functions(cv::ClassCreator<pilot> &cc){

  using namespace cvv8;
  typedef ClassCreator<pilot> CC;

  logstream(LOG_INFO) << "== strapping pilot " << std::flush;
  cc("ping", MethodToInCa<pilot, void (), &pilot::ping>::Call)
    ("destroy", CC::DestroyObjectCallback)
    ("loadGraph",
      MethodToInCa<pilot, void (const std::string&, const std::string&),
        &pilot::load_graph>::Call)
    ("loadSyntheticPowerlaw",
      MethodToInCa<pilot, void (size_t),
        &pilot::load_synthetic_powerlaw>::Call)
    ("transformVertices",
      MethodToInCa<pilot, void (const Handle<Function> &),
        &pilot::transform_vertices>::Call)
    ("fly",
      MethodToInCa<pilot, void (const Handle<Function> &),
        &pilot::fly>::Call)
    ("saveGraph",
      MethodToInCa<
        pilot,
        void (const std::string &, const Handle<Function> &, const Handle<Function> &, bool, bool, bool),
        &pilot::save_graph>::Call);

}

static void bind_global_const (Handle<Object> const & dest){
  using namespace cvv8;
  dest->Set(JSTR("NO_EDGES"), CastToJS(NO_EDGES));
  dest->Set(JSTR("IN_EDGES"), CastToJS(IN_EDGES));
  dest->Set(JSTR("OUT_EDGES"), CastToJS(OUT_EDGES));
  dest->Set(JSTR("ALL_EDGES"), CastToJS(ALL_EDGES));
}

namespace cvv8 {

  CVV8_TypeName_IMPL((pilot), "pilot");

  pilot *ClassCreator_Factory<pilot>::
  Create(Persistent<Object> & jsSelf, Arguments const & argv){
    typedef CtorArityDispatcher<pilotCtors> Proxy;
    pilot * b = Proxy::Call( argv );
    if(b) BMap::Insert( jsSelf, b );
    return b;
  }

  void ClassCreator_Factory<pilot>::
  Delete(pilot *obj){
    BMap::Remove(obj);
    delete obj;
  }

  template <>
  struct ClassCreator_SetupBindings<pilot> {

    static void Initialize(Handle<Object> const & dest) {

      logstream(LOG_INFO) << "== Preparing cockpit " << std::flush;

      // bootstrap class-wrapping code
      typedef ClassCreator<pilot> CC;
      CC & cc (CC::Instance());
      if (cc.IsSealed()) {
        cc.AddClassTo(TypeName<pilot>::Value, dest);
        logstream(LOG_INFO) << "== cockpit ready ==>" << std::endl;
        return;
      }

      bind_member_functions(cc);
      cc.AddClassTo(TypeName<pilot>::Value, dest);
      bind_global_const(dest);

    }

  };

};
