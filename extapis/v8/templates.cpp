#include <graphlab.hpp>
#include "cvv8.hpp"
#include "templates.hpp"
#include "pilot.hpp"

using namespace v8;
using namespace graphlab;

namespace cv = cvv8;

namespace cvv8 {
  typedef pilot::graph_type::vertex_type vertex_type;
  typedef pilot::graph_type::vertex_type edge_type;
};
 
templates::templates(){
  HandleScope handle_scope;
  expose_vertex_type();
  expose_edge_type();
}
    
void templates::expose_vertex_type(){ 
  vertex_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  vertex_templ->SetInternalFieldCount(1);
  vertex_templ->SetAccessor(JSTR("data"), templates::get_vertex_data);
}

void templates::expose_edge_type(){
  edge_templ = Persistent<ObjectTemplate>::New(ObjectTemplate::New());
  edge_templ->SetInternalFieldCount(1);
}

Handle<Value> templates::get_vertex_data(Local<String> property, const AccessorInfo &info) {
  Local<Object> self = info.Holder();
  Local<External> wrap = Local<External>::Cast(self->GetInternalField(0));
  void* ptr = wrap->Value();
  pilot::vertex_data_type value = static_cast<cv::vertex_type*>(ptr)->data();
  return cv::CastToJS(value);
}
