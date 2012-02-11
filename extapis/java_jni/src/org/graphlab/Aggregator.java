package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Aggregator
 * 
 * <p>
 * Aggregates values over all vertices at specified intervals.
 * </p>
 * 
 * @param <V> Vertex type that will be used in {@link #update(Context, Vertex)} 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Aggregator<V extends Vertex> implements Cloneable {

  protected abstract void exec(Context context, V vertex);
  protected abstract void add(Aggregator<V> aggregator);
  protected abstract void finalize(Context context);
  
  /*
   * Required because the aggregator might be executed in parallel.
   * (non-Javadoc)
   * @see java.lang.Object#clone()
   */
  @Override
  protected abstract Aggregator<V> clone();
  
  /**
   * Invoked by proxy aggregator. Creates a pseudo-context and hands it
   * to #exec.
   * @param contextPtr
   * @param vertex
   */
  @SuppressWarnings("unused")
  private void exec(long contextPtr, V vertex){
    if (null == vertex)
      throw new NullPointerException("vertex must not be null."); 
    Context context = new Context(contextPtr);
    exec(context, vertex);
  }
  
  @SuppressWarnings("unused")
  private void finalize(long contextPtr){
    Context context = new Context(contextPtr);
    finalize(context);
  }
  
}
