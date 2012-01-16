package org.graphlab;

import java.util.Map;

/**
 * Updater
 * 
 * <p>The GraphLab engine will invoke an updater on each scheduled node. Extend
 * this class to provide an update function for that node. Note that the update
 * function may update node data, modify edge data, and schedule neighbors, but
 * may not modify the graph structure. You may reuse the updater object on
 * across multiple vertices (this is encouraged).
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Updater {
  
  private final Map<Integer, Integer> mIdMap;
  
  /**
   * Subclasses must call super
   * @param core    core that this updater will be used on
   * @throws IllegalArgumentException
   *          if core is null
   * @throws IllegalStateException
   *          if this updater is constructed before
   *          {@link Core#setGraph(org.graphlab.data.Graph)} is called.      
   */
  public Updater (Core<?> core){
    
    if (null == core) throw new IllegalArgumentException("core must not be null.");
    
    mIdMap = core.idMap();
    if (null == mIdMap)
      throw new IllegalStateException("Must call Core#setGraph before creating updaters.");
    
  }

	/**
	 * Updates the vertex identified by <tt>vertex_id</tt>. Subclasses may wish
	 * to maintain a reference to the graph object.
	 * 
	 * @param context      graphlab context; use {@link Context#schedule(int, Updater)}
	 *                     to schedule vertices.
	 * @param vertexId     application vertex ID
	 */
	public abstract void update(Context context, int vertexId);

  /**
   * When multiple update functors are scheduled to be run on the same function
   * they are added. The default behavior is to simply ignore the later update
   * functors. Override this method to implement your own behavior.
   * @param updater
   */
	public void add(Updater updater){
	  return;
	}
	
	/**
   * Executes the updater on the specified vertex. This is <em>only</em>
   * invoked by the proxy updater in the JNI library.
   * 
   * @param contextPtr
   *        address of graphlab::icontext_type object
   * @param vertexId
   *        application vertex ID
   */
  @SuppressWarnings("unused")
  private void execUpdate (long contextPtr, int vertexId){
    Context context = new Context(contextPtr, mIdMap);
    update(context, vertexId);
  }

}
