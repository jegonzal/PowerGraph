package org.graphlab;

import java.util.Map;
import java.util.NoSuchElementException;

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
  
  /** Reference to the core object that instantiated this context */
  private Core<?> mCore;
  
  /** Address of the associated <tt>graphlab::core</tt> object */
  private long mCorePtr;
  
  /** Address of the associated <tt>graphlab::icontext_type</tt> object */
  private long mContextPtr;
  
  /** Map from application vertex IDs to graphlab vertex IDs */
  private Map<Integer, Integer> mIdMap;

  /**
   * Creates a new context.
   * 
   * @param core
   *          core object that instantiated this context
   * @param corePtr
   *          address of the <tt>graphlab::core</tt> object (must be from
   *          <tt>core</tt>).
   * @param contextPtr
   *          address of the <tt>graphlab::icontext_type</tt> object associated
   *          with this context.
   * @param idMap
   *          map from application vertex IDs to GraphLab vertex IDs
   * @throws NullPointerException
   *           if <tt>core</tt> or <tt>idMap</tt> was null
   */
  protected Context(Core<?> core,
                    long corePtr, long contextPtr,
                    Map<Integer, Integer> idMap) {

    if (null == core || null == idMap)
      throw new NullPointerException("core and idMap must not be null.");

    mCore = core;
    mCorePtr = corePtr;
    mContextPtr = contextPtr;
    mIdMap = idMap;

  }
  
  /**
   * Schedules an update on the specified vertex.
   * 
   * @param vertexId
   *          application vertex ID of vertex to update
   * @param updater
   *          ID of updater to apply on the vertex
   * @throws NullPointerException
   *           if <tt>updater</tt> was null.
   * @throws NoSuchElementException
   *           if <tt>vertexId</tt> did not exist in the graph that was passed
   *           to {@link Core#setGraph(org.graphlab.data.Graph)}.
   */
  public void schedule (int vertexId, Updater updater){
    
    if (null == updater)
      throw new NullPointerException("updater must not be null.");

    Integer glVertexId = mIdMap.get(vertexId);
    if (null == glVertexId)
      throw new NoSuchElementException(
          "vertex did not exist in the graph that was passed to Core#setGraph.");

    // adds updater to core, which creates an ID for the updater
    mCore.addUpdater(updater);
    schedule(mCorePtr, mContextPtr, glVertexId, updater.id());
    
  }
  
  /**
   * Schedules an update on the specified vertex. This invokes the
   * <tt>schedule</tt> method on the associated <tt>graphlab::icontext_type</tt>
   * object (which is accessed using <tt>context_ptr</tt>).
   * 
   * @param core_ptr
   *          address of the associated <tt>graphlab::core</tt> object.
   * @param context_ptr
   *          address of the associated <tt>graphlab::icontext_type</tt> object.
   * @param vertex_id
   *          graphlab vertex ID of vertex to update
   * @param updater_id
   *          ID of updater to apply on the vertex. Updater must be passed to
   *          {@link Core#addUpdater(Updater)} at least once before this method
   *          is invoked.
   */
  private native void schedule(long core_ptr, long context_ptr,
                              int vertex_id, int updater_id);
  
}
