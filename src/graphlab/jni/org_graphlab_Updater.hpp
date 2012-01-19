/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

/**
 * @file org_graphlab_Updater.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#ifndef ORG_GRAPHLAB_UPDATER_HPP
#define ORG_GRAPHLAB_UPDATER_HPP

#include <graphlab.hpp>
#include "org_graphlab_Core.hpp"
#include "org_graphlab_Updater.h"

namespace graphlab {

  /** Proxy edge */
  struct proxy_edge {};
  
  /** Proxy vertex */
  struct proxy_vertex {
    /** corresponding application vertex ID */
    int app_id;
  };
  
  /** Proxy graph */
  typedef graph<proxy_vertex, proxy_edge> proxy_graph;
  
  /**
   * Proxy updater.
   * Mirrors and forwards update calls to the corresponding Java updater.
   * The constructor creates a new reference to the Java object (so that it
   * doesn't get garbage collected.) The destructor will delete the reference
   * to allow the corresponding Java object to be garbaged collected. Note
   * that multiple proxy_updaters may correspond to the same
   * org.graphlab.Updater object.
   */
  class proxy_updater : 
    public iupdate_functor<proxy_graph, proxy_updater> {
    
  private:
    
    /** org.graphlab.Updater object */
    jobject mjava_updater;
    
  public:
  
    /** jni_core type that uses the proxy graph and the proxy updater */
    typedef jni_core<proxy_graph, proxy_updater> core;

    /** context type that uses the proxy graph and the proxy updater */
    typedef iupdate_functor<proxy_graph, proxy_updater>::icontext_type context;
  
    /**
     * Method ID of org.graphlab.Updater#execUpdate.
     */
    static jmethodID java_exec_update;
    
    /**
     * Method ID of org.graphlab.Updater#add.
     */
    static jmethodID java_add;
    
    /**
     * Method ID of org.graphlab.Updater#priority.
     */
    static jmethodID java_priority;
    
    /**
     * Constructor for proxy updater.
     * Initializes this object with the associated Java org.graphlab.Updater
     * object.
     * @param[in] env           JNI environment - used to create a new reference
     *                          to javaUpdater.
     * @param[in] java_updater  Java org.graphlab.Updater object. This constructor
     *                          will create a new reference to the object to prevent
     *                          garbage collection.
     */
    proxy_updater(JNIEnv *env, jobject &java_updater);
    
    /** The default constructor does nothing */
    proxy_updater();
    
    /**
     * Copy constructor for proxy_updater.
     * If \c other has a \c mjava_updater, creates a new reference to it.
     */
    proxy_updater(const proxy_updater& other);
    
    /**
     * Copy assignment operator for proxy_updater.
     * If \c other has a \c mjava_updater, creates a new reference to it.
     */
    proxy_updater &operator=(const proxy_updater &other);
    
    /**
     * Deletes the reference to the Java object so that it may be garbage
     * collected.
     */
    ~proxy_updater();
    
    void operator()(icontext_type& context);
    
    void operator+=(const update_functor_type& other) const;
    
    static void init(JNIEnv *env);
    
  };
  
}

#endif