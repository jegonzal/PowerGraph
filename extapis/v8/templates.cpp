#include <graphlab.hpp>
#include "templates.hpp"

using namespace v8;
using namespace graphlab;

namespace cv = cvv8;

templates::templates(){
  HandleScope handle_scope;
  expose_vertex_type();
  expose_edge_type();
  expose_context_type();
}

void templates::expose_vertex_type(){

  vertex_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  vertex_templ->SetInternalFieldCount(2); // see JSToNative<vertex_type>

  // bind vertex.data() as a JS property
  vertex_templ->SetAccessor(JSTR("data"), templates::get_vertex_data, templates::set_vertex_data);

  // bind vertex.num_out_edges as JS member function
  vertex_templ->Set(JSTR("numOutEdges"), FunctionTemplate::New(templates::get_vertex_num_out_edges));

  // bind vertex.id as JS member function
  vertex_templ->Set(JSTR("id"), FunctionTemplate::New(templates::get_vertex_id));

}

void templates::expose_edge_type(){

  edge_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  edge_templ->SetInternalFieldCount(1);

  // bind source() method
  edge_templ->Set(JSTR("source"), FunctionTemplate::New(
    cv::ConstMethodToInCa<
      cv::edge_type,
      cv::vertex_type (),
      &cv::edge_type::source
    >::Call
  ));

  // bind target() method
  edge_templ->Set(JSTR("target"), FunctionTemplate::New(
    cv::ConstMethodToInCa<
      cv::edge_type,
      cv::vertex_type (),
      &cv::edge_type::target
    >::Call
  ));

}

void templates::expose_context_type(){

  context_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  context_templ->SetInternalFieldCount(1);

  // bind signal() method
  context_templ->Set(JSTR("signal"), FunctionTemplate::New(
    cv::MethodToInCaVoid<
      cv::context_type,
      void (
        const cv::context_type::vertex_type&,
        const cv::message_type&
      ),
      &cv::context_type::signal
    >::Call
  ));

}

templates::~templates(){
  // dispose persistent handles
  vertex_templ.Dispose();
  edge_templ.Dispose();
  context_templ.Dispose();
}

Handle<Value> templates::get_vertex_data(Local<String> property, const AccessorInfo &info) {
  HandleScope handle_scope;
  const Handle<Value> &h = info.Holder();
  return handle_scope.Close(cv::CastToJS(cv::JSToNative<cv::vertex_type>()(h).data()));
}

void templates::set_vertex_data(Local<String> property, Local<Value> value, const AccessorInfo& info){
  const Handle<Value> &h = info.Holder();
  cv::JSToNative<cv::vertex_type>()(h).data() = cv::CastFromJS<pilot::vertex_data_type>(value);
}

Handle<Value> templates::get_vertex_num_out_edges(const Arguments &argv){
  HandleScope handle_scope;
  // cannot use CastFromJS here because there vertex_type casting is unconventional
  const cv::vertex_type vertex = cv::JSToNative<cv::vertex_type>()(argv.This());
  return handle_scope.Close(cv::CastToJS(vertex.num_out_edges()));
}

Handle<Value> templates::get_vertex_id(const Arguments &argv){
  HandleScope handle_scope;
  const cv::vertex_type vertex = cv::JSToNative<cv::vertex_type>()(argv.This());
  return handle_scope.Close(cv::CastToJS(vertex.id()));
}

namespace cvv8 {

  CVV8_TypeName_IMPL((vertex_type), "vertex");
  CVV8_TypeName_IMPL((edge_type), "edge");
  CVV8_TypeName_IMPL((context_type), "context");

  Handle<Value> NativeToJS<vertex_type>::operator()(const vertex_type &vertex){
    Local<Object> v = pilot::get_templates().vertex_templ->NewInstance();
    v->SetInternalField(0, External::New((void *) &vertex.graph_ref));
    v->SetInternalField(1, cv::CastToJS(vertex.local_id()));
    return v;
  }

  Handle<Value> NativeToJS<edge_type>::operator()(const edge_type &edge){
    Local<Object> e = pilot::get_templates().edge_templ->NewInstance();
    e->SetInternalField(0, External::New((void *) &edge));
    return e;
  }

  Handle<Value> NativeToJS<context_type>::operator()(const context_type &context){
    Local<Object> c = pilot::get_templates().context_templ->NewInstance();
    c->SetInternalField(0, External::New((void *) &context));
    return c;
  }

  Handle<Value> NativeToJS<edge_dir_type>::operator()(const edge_dir_type &edge_dir){
    return CastToJS(static_cast<int32_t>(edge_dir));
  }

  vertex_type JSToNative<vertex_type>::operator()(const Handle<Value> &h) const {
    const Local<Object> &obj(Object::Cast(*h));
    graph_type *graph = (graph_type *) obj->GetPointerFromInternalField(0);
    graph_type::lvid_type lvid = CastFromJS<graph_type::lvid_type>(obj->GetInternalField(1));
    return cv::vertex_type(*graph, lvid);
  }

  edge_dir_type JSToNative<edge_dir_type>::operator()(const Handle<Value> &h) const {
    return static_cast<edge_dir_type>(CastFromJS<int32_t>(h));
  }

  graphlab::empty JSToNative<graphlab::empty>::operator()(const Handle<Value> &h) const {
    return graphlab::empty();
  }

};
