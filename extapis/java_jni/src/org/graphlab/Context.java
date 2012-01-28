package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * GraphLab Context.
 * 
 * <p>
 * This mirrors <tt>graphlab::icontext_type</tt>. Applications should
 * <em>never</em> instantiate this class; instances of this class will be passed
 * to the {@link org.graphlab.Updater updater} when the updater is invoked by
 * the GraphLab scheduler. The updater may then use the context object to
 * schedule updates on vertices.
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public final class Context {
  
  /** Address of the associated <tt>graphlab::icontext_type</tt> object */
  private long mContextPtr;
  
  /**
   * Creates a new context.
   * @param contextPtr
   *          address of the <tt>graphlab::icontext_type</tt> object associated
   *          with this context.
   * @throws NullPointerException
   *           if <tt>core</tt> or <tt>idMap</tt> was null
   */
  protected Context(long contextPtr) {
    mContextPtr = contextPtr;
  }
  
  /**
   * Schedules an update on the specified vertex.
   * 
   * @param vertex
   *          vertex to update
   * @param updater
   *          ID of updater to apply on the vertex
   * @throws NullPointerException
   *           if <tt>updater</tt> or <tt>vertex</tt> was null.
   */
  public void schedule (Vertex vertex, Updater<?> updater){
    
    if (null == updater || null == vertex)
      throw new NullPointerException("updater must not be null.");

    // adds updater to core, which creates an ID for the updater
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
  private native void schedule(long context_ptr, Updater<?> updater, int vertexId);
  
}
