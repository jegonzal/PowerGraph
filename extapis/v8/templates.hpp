/**
 * @file templates.hpp
 */

#ifndef GRAPHLAB_TEMPLATES_HPP
#define GRAPHLAB_TEMPLATES_HPP

#include <v8.h>

namespace graphlab {
   
  /**
   * Defines object templates for vertex and edge.
   * TODO: figure out a way to do this in cvv8
   */
  class templates {
  public:
    v8::Persistent<v8::ObjectTemplate> vertex_templ;
    v8::Persistent<v8::ObjectTemplate> edge_templ;
  public:
    templates();
    static v8::Handle<v8::Value> get_vertex_data(v8::Local<v8::String> property, const v8::AccessorInfo &info);
  private:
    /** Adds a template for vertex_type */
    void expose_vertex_type();
    /** Adds a template for edge_type */
    void expose_edge_type();
  };

};

#endif // GRAPHLAB_TEMPLATES_HPP
