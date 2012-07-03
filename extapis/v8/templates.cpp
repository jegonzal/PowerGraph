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
  vertex_templ->SetInternalFieldCount(2);
  vertex_templ->SetAccessor(JSTR("data"), templates::get_vertex_data, templates::set_vertex_data);
  vertex_templ->SetAccessor(JSTR("num_out_edges"), templates::get_vertex_num_out_edges);
  vertex_templ->SetAccessor(JSTR("id"), templates::get_vertex_id);
  // TODO add other getters
}

void templates::expose_edge_type(){
  edge_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  edge_templ->SetInternalFieldCount(1);
  edge_templ->SetAccessor(JSTR("source"), templates::get_edge_source);
  edge_templ->SetAccessor(JSTR("target"), templates::get_edge_target);
}

void templates::expose_context_type(){
  context_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  context_templ->SetInternalFieldCount(1);
  context_templ->Set(JSTR("signal"), FunctionTemplate::New(
    cv::MethodToInCaVoid<
      cv::context_type,
      void (
        const cv::context_type::vertex_type&,
        const cv::message_type &message
      ),
      &cv::context_type::signal
    >::Call
  ));
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

Handle<Value> templates::get_vertex_num_out_edges(Local<String> property, const AccessorInfo &info){
  HandleScope handle_scope;
  const Handle<Value> &h = info.Holder();
  const cv::vertex_type &vertex = cv::CastFromJS<const cv::vertex_type &>(h);
  return handle_scope.Close(cv::CastToJS(vertex.num_out_edges()));
}

Handle<Value> templates::get_vertex_id(Local<String> property, const AccessorInfo &info){
  HandleScope handle_scope;
  const Handle<Value> &h = info.Holder();
  const cv::vertex_type &vertex = cv::CastFromJS<const cv::vertex_type &>(h);
  return handle_scope.Close(cv::CastToJS(vertex.id()));
}

// TODO: convert to method invocation
Handle<Value> templates::get_edge_source(Local<String> property, const AccessorInfo &info){
  HandleScope handle_scope;
  Local<Object> self = info.Holder();
  Local<External> wrap = Local<External>::Cast(self->GetInternalField(0));
  void *ptr = wrap->Value();
  cv::vertex_type vertex = static_cast<cv::edge_type*>(ptr)->source();
  return handle_scope.Close(cv::CastToJS(vertex));
}

// TODO: convert to method invocation
Handle<Value> templates::get_edge_target(Local<String> property, const AccessorInfo &info){
  HandleScope handle_scope;
  Local<Object> self = info.Holder();
  Local<External> wrap = Local<External>::Cast(self->GetInternalField(0));
  void *ptr = wrap->Value();
  cv::vertex_type vertex = static_cast<cv::edge_type*>(ptr)->target();
  return handle_scope.Close(cv::CastToJS(vertex));
}

namespace cvv8 {

  // cannot use CVV8TypeName_IMPL because context_type is abstract
  const char *TypeName<context_type>::Value = "context_type";

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

  vertex_type JSToNative<vertex_type>::operator()(const Handle<Value> &h) const {
    const Local<Object> &obj(Object::Cast(*h));
    graph_type *graph = (graph_type *) obj->GetPointerFromInternalField(0);
    graph_type::lvid_type lvid = CastFromJS<graph_type::lvid_type>(obj->GetInternalField(1));
    return cv::vertex_type(*graph, lvid);
  }

  graphlab::empty JSToNative<const graphlab::empty &>::operator()(const Handle<Value> &h) const {
    return graphlab::empty();
  } 

};
