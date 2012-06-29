/**
 * @file templates.hpp
 */

#ifndef GRAPHLAB_TEMPLATES_HPP
#define GRAPHLAB_TEMPLATES_HPP

#include <v8.h>
#include "cvv8.hpp"
#include "pilot.hpp"

namespace graphlab {
   
  /**
   * Defines object templates for vertex_type and edge_type.
   * TODO: refactor to use cvv8 and CC.
   */
  class templates {
  public:
    /** Object template for vertex_type */
    v8::Persistent<v8::ObjectTemplate> vertex_templ;
    /** Object template for edge_type */
    v8::Persistent<v8::ObjectTemplate> edge_templ;
    /** Object template for context_type */
    v8::Persistent<v8::ObjectTemplate> context_templ;
  public:
    templates();
    static v8::Handle<v8::Value> get_vertex_data(v8::Local<v8::String> property, const v8::AccessorInfo &info);
    static void set_vertex_data(v8::Local<v8::String> property, v8::Local<v8::Value> value, const v8::AccessorInfo& info);
    static v8::Handle<v8::Value> get_vertex_num_out_edges(v8::Local<v8::String> property, const v8::AccessorInfo &info);
    static v8::Handle<v8::Value> get_vertex_id(v8::Local<v8::String> property, const v8::AccessorInfo &info);
    static v8::Handle<v8::Value> get_edge_source(v8::Local<v8::String> property, const v8::AccessorInfo &info);
    static v8::Handle<v8::Value> get_edge_target(v8::Local<v8::String> property, const v8::AccessorInfo &info);
  private:
    /** Adds an object template for vertex_type */
    void expose_vertex_type();
    /** Adds an object template for edge_type */
    void expose_edge_type();
    /** Adds an object template for context type */
    void expose_context_type();
  };

};

namespace cvv8 {

  typedef graphlab::pilot::graph_type graph_type;
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type edge_type;
  typedef graphlab::js_proxy::icontext_type context_type;
  typedef graphlab::js_proxy::message_type message_type;

  CVV8_TypeName_DECL((context_type));

  /**
   * Convenience functor for casting from vertex to a JS wrapper.
   * Saves pointer to graph and vertex id in the JS object.
   */
  template<>
  struct NativeToJS<vertex_type> {
    v8::Handle<v8::Value> operator()(const vertex_type &vertex);
  };

  /**
   * Convenience functor for casting from edge to a JS wrapper.
   * Note that the edge must exist for the lifetime of the returned
   * handle.
   */
  template<>
  struct NativeToJS<edge_type> {
    v8::Handle<v8::Value> operator()(const edge_type &edge);
  };

  /**
   * Convenience functor for casting from context to JS wrapper.
   * Note that the context must exist for the lifetime of the returned
   * handle.
   */
  template <>
  struct NativeToJS<context_type> {
    v8::Handle<v8::Value> operator()(const context_type &context);
  };

  template <>
  struct JSToNative<vertex_type> {
    vertex_type operator()(const v8::Handle<v8::Value> &h) const;
  };

  /** Not sure why I need this but I do */
  template <>
  struct JSToNative<const vertex_type &> :
    JSToNative<vertex_type>{
    typedef vertex_type ResultType;
  };
  
  /** Convertes js object to context using pointer stored in
   * its internal field. */
  template <>
  struct JSToNative<context_type> :
    JSToNative_ObjectWithInternalFields<context_type>{};

  /** Convertes js object to graph using pointer stored in
   * its internal field. */
  template <>
  struct JSToNative<graph_type> :
    JSToNative_ObjectWithInternalFields<graph_type>{};

  /** Convertes js object to graph using pointer stored in
   * its internal field. TODO: this is broken */
  template <>
  struct JSToNative<message_type> :
    JSToNative_ObjectWithInternalFields<message_type>{};

  /** Special case for graphlab::empty */
  template <>
  struct JSToNative<graphlab::empty const &> {
    typedef graphlab::empty ResultType;
    ResultType operator()(const v8::Handle<v8::Value> &h) const;
  };

};

#endif // GRAPHLAB_TEMPLATES_HPP
