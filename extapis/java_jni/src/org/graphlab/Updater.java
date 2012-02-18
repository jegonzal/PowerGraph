package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Updater
 * 
 * <p>
 * The GraphLab engine will invoke an updater on each scheduled node. Extend
 * this class to provide an update function for that node. Note that the update
 * function may update node data, modify edge data, and schedule neighbors, but
 * may not modify the graph structure. You may reuse the updater object on
 * across multiple vertices (this is encouraged).
 * </p>
 * 
 * @param <V> Vertex type that will be used in {@link #update(Context, Vertex)} 
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Updater<V extends Vertex> {

  /**
   * Updates the vertex. Subclasses may wish to
   * maintain a reference to the graph object.
   * 
   * @param context
   *          graphlab context; use {@link Context#schedule(Vertex, Updater)} to
   *          schedule vertices.
   * @param vertex
   *          vertex to be updated
   */
  protected abstract void update(Context context, V vertex);

  /**
   * When multiple update functors are scheduled to be run on the same function
   * they are added. The default behavior is to simply ignore the later update
   * functors. Override this method to implement your own behavior.
   * 
   * @param updater
   */
  protected void add(Updater<V> updater) {
    return;
  }
  
  /**
   * Get the priority of the update functor. Defaults to 0.
   * @return priority
   */
  public double priority(){
    return 0;
  }
  
  /*
   * Required because multiple updaters might be executed in parallel.
   * (non-Javadoc)
   * @see java.lang.Object#clone()
   */
  @Override
  protected abstract Updater<V> clone();

  /**
   * Executes the updater on the specified vertex. This is <em>only</em> invoked
   * by the proxy updater in the JNI library.
   * 
   * @param contextPtr
   *          address of graphlab::icontext_type object
   * @param vertexId
   *          application vertex ID
   */
  @SuppressWarnings("unused")
  private void update(long contextPtr, V vertex) {
    
    // for debugging purposes, throw up on bogus vertices from GraphLab
    if (null == vertex)
       throw new NullPointerException("vertex must not be null.");
      
    Context context = new Context(contextPtr);
    update(context, vertex);
    
  }

}
