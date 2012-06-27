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
  vertex_templ->SetAccessor(JSTR("data"), templates::get_vertex_data);
  vertex_templ->SetAccessor(JSTR("num_out_edges"), templates::get_vertex_num_out_edges);
  // TODO add other getters
}

void templates::expose_edge_type(){
  edge_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  edge_templ->SetInternalFieldCount(1);
  edge_templ->SetAccessor(JSTR("source"), templates::get_edge_source);
}

void templates::expose_context_type(){
  context_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  context_templ->SetInternalFieldCount(1);
  context_templ->Set(JSTR("signal"), FunctionTemplate::New(
    cv::MethodToInCa<
      cv::context_type,
      void (const cv::context_type::vertex_type&, const cv::context_type::message_type&),
      &cv::context_type::signal
    >::Call
  ));
}

pilot::graph_type::vertex_type templates::castToVertex(const v8::AccessorInfo &info){
  Local<Object> self = info.Holder();
  Local<External> graph_wrap = Local<External>::Cast(self->GetInternalField(0));
  pilot::graph_type *graph_ref = (pilot::graph_type *) graph_wrap->Value();
  pilot::graph_type::vertex_id_type id = cv::CastFromJS<pilot::graph_type::vertex_id_type>(self->GetInternalField(1));
  // rebuild vertex from graph reference and ID
  return cv::vertex_type(*graph_ref, id);
}

Handle<Value> templates::get_vertex_data(Local<String> property, const AccessorInfo &info) {
  return cv::CastToJS(castToVertex(info).data());
}

Handle<Value> templates::get_vertex_num_out_edges(Local<String> property, const AccessorInfo &info){
  return cv::CastToJS(castToVertex(info).num_out_edges());
}

Handle<Value> templates::get_edge_source(Local<String> property, const AccessorInfo &info){
  Local<Object> self = info.Holder();
  Local<External> wrap = Local<External>::Cast(self->GetInternalField(0));
  void *ptr = wrap->Value();
  cv::vertex_type vertex = static_cast<cv::edge_type*>(ptr)->source();
  return cv::CastToJS(vertex);
}

namespace cvv8 {

  // cannot use CVV8TypeName_IMPL because context_type is abstract
  const char *TypeName<context_type>::Value = "context_type";

  Handle<Value> NativeToJS<vertex_type>::operator()(const vertex_type vertex){
    Local<Object> v = pilot::get_templates().vertex_templ->NewInstance();
    v->SetInternalField(0, External::New((void *) &vertex.graph_ref));
    v->SetInternalField(1, cv::CastToJS(vertex.id()));
    return v;
  }
  
  Handle<Value> NativeToJS<vertex_type *>::operator()(const vertex_type *vertex){
    Local<Object> v = pilot::get_templates().vertex_templ->NewInstance();
    v->SetInternalField(0, External::New((void *) vertex));
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
    graph_type::vertex_id_type id = CastFromJS<graph_type::vertex_id_type>(obj->GetInternalField(1));
    return cv::vertex_type(*graph, id);
  }

};
