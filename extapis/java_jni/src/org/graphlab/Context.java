package org.graphlab;

import java.util.Map;
import java.util.NoSuchElementException;

public final class Context {
  
  /**
   * Thread-safe singleton container.
   * This class is loaded only when {@link Context#getInstance} is first
   * accessed and is loaded exactly once (Java language guarantee).
   * @see <a href="http://www.cs.umd.edu/~pugh/java/memoryModel/DoubleCheckedLocking.html">Singleton</a>
   */
  private static class SingletonHolder {
    public static final Context singleton = new Context();
  }
  
  private Map<Integer, Integer> mIdMap = null;
  
  private Context (){
  }
  
  public static Context getInstance (){
    return SingletonHolder.singleton;
  }
  
  protected void setIdMap (Map<Integer, Integer> map){
    mIdMap = map;
  }
  
  public void schedule (long corePtr, long contextPtr, int vertexId, Updater updater){
    
    if (null == updater)
      throw new NullPointerException("updater must not be null.");
    
    if (null == mIdMap)
      throw new NullPointerException("mIdMap not initialized.");

    Integer glVertexId = mIdMap.get(vertexId);
    if (null == glVertexId)
      throw new NoSuchElementException("vertex did not exist in the graph that was passed to Core#setGraph.");

    schedule(corePtr, contextPtr, glVertexId, updater.id());
    
  }
  
  private native void schedule (long core_ptr, long context_ptr, int vertex_id, int updater_id);
  
}
