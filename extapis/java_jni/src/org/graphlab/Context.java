package org.graphlab;

import java.util.Map;
import java.util.NoSuchElementException;

public final class Context {
  
  private Core<?> mCore;
  private long mCorePtr;
  private long mContextPtr;
  
  /** Map from application vertex IDs to graphlab vertex IDs */
  private Map<Integer, Integer> mIdMap;
  
  protected Context (Core<?> core, long corePtr, long contextPtr, Map<Integer, Integer> idMap){
    
    if (null == core || null == idMap)
      throw new NullPointerException ("core and idMap must not be null.");
    
    mCore = core;
    mCorePtr = corePtr;
    mContextPtr = contextPtr;
    mIdMap = idMap;
     
  }
  
  /**
   * Schedules the specified vertex for updating.
   * @param vertexId
   *        application vertex ID of vertex to update
   * @param updater
   *        ID of updater to apply on the vertex
   */
  public void schedule (int vertexId, Updater updater){
    
    if (null == updater)
      throw new NullPointerException("updater must not be null.");

    Integer glVertexId = mIdMap.get(vertexId);
    if (null == glVertexId)
      throw new NoSuchElementException("vertex did not exist in the graph that was passed to Core#setGraph.");

    mCore.addUpdater(updater);
    schedule(mCorePtr, mContextPtr, glVertexId, updater.id());
    
  }
  
  /**
   * Schedules the specified vertex for updating.
   * @param core_ptr
   * @param context_ptr
   * @param vertex_id
   *        graphlab vertex ID of vertex to update
   * @param updater_id
   *        ID of updater to apply on the vertex
   */
  private native void schedule (long core_ptr, long context_ptr, int vertex_id, int updater_id);
  
}
