package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Aggregator
 * 
 * <p>
 * Aggregates values over all vertices at specified intervals. Analogous to performing
 * map-reduce on the graph.
 * </p>
 * 
 * @param <V> vertex type that will be used in {@link #exec(Context, Vertex)}.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Aggregator<V extends Vertex, A extends Aggregator<V, A>> implements Cloneable {

  static { initNative(); }
  
  /**
   * Executes operation on a single vertex.
   * @param context
   * @param vertex
   */
  protected abstract void exec(Context context, V vertex);
  
  /**
   * Merges results of multiple aggregators.
   * @param aggregator
   */
  protected abstract void add(A aggregator);
  
  /**
   * Called when aggregation is completed. Usually used to output results.
   * @param context
   */
  protected abstract void finalize(Context context);
  
  /*
   * Required because the aggregator might be executed in parallel.
   * (non-Javadoc)
   * @see java.lang.Object#clone()
   */
  @Override
  protected abstract A clone();
  
  /**
   * Invoked by proxy aggregator. Creates a pseudo-context and hands it
   * to {@link #exec(Context, Vertex)}.
   * @param contextPtr
   * @param vertex
   */
  private void exec(long contextPtr, V vertex){
    if (null == vertex)
      throw new NullPointerException("vertex must not be null."); 
    Context context = new Context(contextPtr);
    exec(context, vertex);
  }
  
  /**
   * Invoked by the proxy aggregator. Creates a pseudo-context and hands it
   * to {@link #finalize()}.
   * @param contextPtr
   */
  private void finalize(long contextPtr){
    Context context = new Context(contextPtr);
    finalize(context);
  }
  
  /**
   * Initialize native class (set field IDs and method IDs)
   */
  private static native void initNative();
  
}
