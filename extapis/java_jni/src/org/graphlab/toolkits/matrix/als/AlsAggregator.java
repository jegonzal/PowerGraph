package org.graphlab.toolkits.matrix.als;

import org.apache.log4j.Logger;
import org.graphlab.Aggregator;
import org.graphlab.Context;

/**
 * Aggregate RMSE of every vertex.
 * 
 * @author Jiunn Haur Lim
 */
public class AlsAggregator
extends Aggregator<AlsVertex, AlsAggregator> {
  
  private static final Logger logger = Logger.getLogger(AlsAggregator.class);

  private double mSSE = 0;
  private AlsGraph mGraph;

  protected AlsAggregator(AlsGraph graph) {
    mGraph = graph;
  }

  @Override
  protected void exec(Context context, AlsVertex vertex) {
    int numEdges = mGraph.edgesOf(vertex).size();
    if (0 == numEdges) return;
    mSSE += vertex.mSSE;
  }

  @Override
  protected void add(AlsAggregator other) {
    mSSE += other.mSSE;
  }

  @Override
  protected void finalize(Context context) {
    logger.info("-------- Train Results --------");
    logger.info("Average RMS: " + Math.sqrt(mSSE / (2 * mGraph.edgeSet().size())));
    mSSE = 0;
  }

  @Override
  protected AlsAggregator clone() {
    AlsAggregator agg = new AlsAggregator(mGraph);
    agg.mSSE = mSSE;
    return agg;
  }

}
