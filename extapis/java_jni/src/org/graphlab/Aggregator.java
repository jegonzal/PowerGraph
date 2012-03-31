package org.graphlab;

import org.graphlab.data.Vertex;

/**
 * Aggregates values over all vertices at specified intervals. This is analogous
 * to performing map-reduce on the graph.
 * 
 * <p>
 * To use, provide an implementation that overrides
 * {@link #exec(Context, Vertex)}, {@link #add(Aggregator)}, {@link #clone()},
 * and {@link #finalize(long)}. Then, pass an instance of your aggregator to
 * {@link Core#addAggregator(String, Aggregator, long)}, specifying the
 * frequency at which you want the aggregator to be initiated. When your
 * aggregator is run, <tt>exec</tt> will be invoked on every vertex, and
 * <tt>add</tt> will be invoked to merge the results of two aggregators.
 * </p>
 * 
 * <h3>Generics</h3>
 * <p>
 * Most of the time, <tt>V</tt> will take the type of vertex that your graph has
 * and <tt>A</tt> will take the type of your aggregator. Your class signature
 * should look like the following:
 * </p>
 * <pre>
 *   private static class Agg extends Aggregator&lt;AlsVertex, Agg&gt;
 * </pre>
 * <p>We recommend that you follow this pattern closely, unless you are very
 * familiar with generics.</p>
 * 
 * @param <V>
 *          Vertex type that will be used in {@link #exec(Context, Vertex)}.
 * @param <A>
 *          For self-templating.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public abstract class Aggregator<V extends Vertex, A extends Aggregator<V, A>>
    implements Cloneable {

  // initialize method IDs
  static { initNative(); }
  
  /**
   * Executes operation on a single vertex.
   * @param context
   * @param vertex
   *          The vertex from which data may be collected.
   *          
   */
  protected abstract void exec(Context context, V vertex);
  
  /**
   * Merges results of multiple aggregators.
   * @param aggregator
   *          Another aggregator whose results may be merged with
   *          those of this aggregator.
   */
  protected abstract void add(A aggregator);
  
  /**
   * Called when aggregation is completed. Usually used to output results.
   * @param context
   */
  protected abstract void finalize(Context context);
  
  /**
   * Clones the aggregator.
   * 
   * Aggregation is often executed in parallel by making copies of existing
   * aggregators.
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
