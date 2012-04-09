package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Execution context.
 * 
 * <p>Applications should <em>never</em> instantiate this class; instances of this
 * class will be passed to the {@link org.graphlab.Updater updater} when the
 * updater is invoked by the GraphLab scheduler. The updater may then use the
 * context object to schedule updates on vertices.</p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public final class Context {
  
  /** Address of the associated <tt>graphlab::icontext_type</tt> object */
  private long mContextPtr;
  
  /**
   * Creates a new context.
   * 
   * <p>Applications should <em>never</em> instantiate this class.</p>
   * @param contextPtr
   *          address of the <tt>graphlab::icontext_type</tt> object associated
   *          with this context.
   */
  protected Context(long contextPtr) {
    mContextPtr = contextPtr;
  }
  
  /**
   * Schedules an update on the given vertex.
   * 
   * <p>When GraphLab invokes the updater on a vertex, it passes the
   * updater the context (this object) and the vertex. When the updater
   * has completed its computation, it may use the context object to
   * schedule updates on neighboring vertices (by calling this method.)</p>
   * 
   * @param vertex
   *          vertex to update
   * @param updater
   *          ID of updater to apply on the vertex
   * @throws NullPointerException
   *           if <tt>updater</tt> or <tt>vertex</tt> was null.
   */
  public void schedule (Vertex vertex, Updater<?, ?, ?> updater){
    
    if (null == updater || null == vertex)
      throw new NullPointerException("updater and vertex must not be null.");

    // schedules an update on the vertex using the proxy ID
    schedule(mContextPtr, updater, vertex.rawId());
    
  }
  
  /**
   * Schedules an update on the specified vertex. This invokes the
   * <tt>schedule</tt> method on the associated <tt>graphlab::icontext_type</tt>
   * object (which is accessed using <tt>context_ptr</tt>).
   * 
   * @param context_ptr
   *          address of the associated <tt>graphlab::icontext_type</tt> object.
   * @param updater
   * @param vertex
   *          graphlab id of vertex to update
   */
  private native void schedule(long context_ptr, Updater<?, ?, ?> updater, int vertexId);
  
}
