package org.graphlab.toolkits.matrix.als;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

/**
 * Graph Type for ALS.
 * @author Jiunn Haur Lim
 */
public class AlsGraph extends SimpleWeightedGraph<AlsVertex, DefaultWeightedEdge> {

  private static final long serialVersionUID = -5986194133580275107L;

  public AlsGraph(Class<? extends DefaultWeightedEdge> edgeClass) {
    super(edgeClass);
  }

}
