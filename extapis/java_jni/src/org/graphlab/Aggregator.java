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
 * <h3>Example</h3>
 * <p>The complete source for this example is available as {@link org.graphlab.demo.Sum}.</p>
 * <p>Suppose we would like to sum over all vertices in the graph. To do so, we
 * implement an aggregator as follows:</p>
<pre>
  private static class SumAggregator extends Aggregator<ScalarVertex, SumAggregator> {

    private int sum;
    
    protected void add(SumAggregator agg) {
      sum += agg.sum;
    }

    protected SumAggregator clone() {
      SumAggregator agg = new SumAggregator();
      agg.sum = sum;
      return agg;
    }

    protected void exec(Context context, ScalarVertex vertex) {
      sum += vertex.value();
    }

    protected void finalize(Context context) {
      System.out.println("Sum: " + sum);
    }
    
  }
</pre>
 * <p>Then, to add it to the core and manually execute it, we use</p>
<pre>
  core.addAggregator("agg", new SumAggregator(), 0);
  core.aggregateNow("agg");
</pre>
 * <p>Instead of using 0 in the third parameter, you may also specify the
 * frequency at which the aggregator should be executed. For more information,
 * refer to {@link Core#addAggregator}.</p>
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
   * 
   * When the core performs an aggregation, the following happens:
   * <ol>
   *  <li>The associated aggregator is cloned (number of clones may vary.)</li>
   *  <li><code>exec</code> is invoked once per vertex.</li>
   *  <li>{@link #add(Aggregator)} is invoked to merge the results of pairs of aggregators.</li>
   *  <li>When there are no more vertices to <code>exec</code> and no pairs of aggregators
   *  to merge, {@link #finalize(Context)} is invoked.</li>
   * </ol>
   * 
   * @param context
   * @param vertex
   *          the vertex from which data may be collected.      
   */
  protected abstract void exec(Context context, V vertex);
  
  /**
   * Merges results of multiple aggregators.
   * 
   * When the core performs an aggregation, the following happens:
   * <ol>
   *  <li>The associated aggregator is cloned (number of clones may vary.)</li>
   *  <li>{@link #exec(Context, Vertex)} is invoked once per vertex.</li>
   *  <li><code>add</code> is invoked to merge the results of pairs of aggregators.</li>
   *  <li>When there are no more vertices to <code>exec</code> and no pairs of aggregators
   *  to merge, {@link #finalize(Context)} is invoked.</li>
   * </ol>
   * 
   * @param aggregator
   *          another aggregator whose results may be merged with
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
   * @see #add(Aggregator)
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
