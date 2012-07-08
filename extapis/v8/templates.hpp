/**
 * @file templates.hpp
 * Javascript bindings for vertex_type, edge_type, message_type, and context_type.
 */

#ifndef GRAPHLAB_TEMPLATES_HPP
#define GRAPHLAB_TEMPLATES_HPP

#include <v8.h>
#include "cvv8.hpp"
#include "pilot.hpp"

namespace graphlab {

  /**
   * Defines object templates for vertex_type, edge_type, and context_type, and context_type.
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
    virtual ~templates();
    /** Pretend that vertex.data() is a property - getter */
    static v8::Handle<v8::Value> get_vertex_data(v8::Local<v8::String> property, const v8::AccessorInfo &info);
    /** Pretend that vertex.data() is a property - setter */
    static void set_vertex_data(v8::Local<v8::String> property, v8::Local<v8::Value> value, const v8::AccessorInfo& info);
    /** vertex.num_out_edges */
    static v8::Handle<v8::Value> get_vertex_num_out_edges(const v8::Arguments &argv);
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

  CVV8_TypeName_DECL((vertex_type));
  CVV8_TypeName_DECL((edge_type));
  CVV8_TypeName_DECL((context_type));

  /**
   * Convenience functor for casting from vertex to a JS wrapper.
   * Saves pointer to graph and vertex id in the JS object.
   * @internal
   */
  template<>
  struct NativeToJS<vertex_type> {
    v8::Handle<v8::Value> operator()(const vertex_type &vertex);
  };

  /**
   * Convenience functor for casting from edge to a JS wrapper.
   * Note that the edge must exist for the lifetime of the returned
   * handle.
   * @internal
   */
  template<>
  struct NativeToJS<edge_type> {
    v8::Handle<v8::Value> operator()(const edge_type &edge);
  };

  /**
   * Convenience functor for casting from context to JS wrapper.
   * Note that the context must exist for the lifetime of the returned
   * handle.
   * @internal
   */
  template <>
  struct NativeToJS<context_type> {
    v8::Handle<v8::Value> operator()(const context_type &context);
  };

  template <>
  struct JSToNative<vertex_type> {
    vertex_type operator()(const v8::Handle<v8::Value> &h) const;
  };

  /**
   * The default implementation casts the JS type to a pointer and uses that
   * as a reference to the vertex object. That doesn't work for us because:
   * - edge.source() returns by value, which is stored on the stack. I cannot
   *   return this to the JS layer because the object might be used somewhere
   *   else. Therefore, NativeToJS creates a new JS instance that holds the
   *   graph ref and the local ID.
   * - edge.signal() takes a vertex reference, which uses this functor to
   *   cast it to a C++ object. The default implementation is incompatible w.
   *   what I do in NativeToJS<vertex_type>.
   * @internal
   */
  template <>
  struct JSToNative<const vertex_type &> :
    JSToNative<vertex_type> {
    typedef vertex_type ResultType;
  };

  /**
   * Converts js object to context using pointer stored in its internal field.
   * @internal
   */
  template <>
  struct JSToNative<context_type> :
    JSToNative_ObjectWithInternalFields<context_type>{};

  /**
   * Converts js object to edge using pointer stored in its internal field.
   * @internal
   */
  template <>
  struct JSToNative<edge_type> :
    JSToNative_ObjectWithInternalFields<edge_type>{};

  /**
   * Converts js object to graph using pointer stored in its internal field.
   * @internal
   */
  template <>
  struct JSToNative<graph_type> :
    JSToNative_ObjectWithInternalFields<graph_type>{};

  // TODO: handle more generic messages

  /** Special case for graphlab::empty */
  template <>
  struct JSToNative<graphlab::empty> {
    graphlab::empty operator()(const v8::Handle<v8::Value> &h) const;
  };

  /** Special case for graphlab::empty */
  template <>
  struct JSToNative<const graphlab::empty &> :
    JSToNative<graphlab::empty> {
    typedef graphlab::empty ResultType;
  };

};

#endif // GRAPHLAB_TEMPLATES_HPP
