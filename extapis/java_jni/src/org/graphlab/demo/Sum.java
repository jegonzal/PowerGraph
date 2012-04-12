package org.graphlab.demo;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Aggregator;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Sums over all vertices.
 * 
 * Before running, add "-Djava.library.path=${resource_loc:/GraphLab-App/lib}" to the
 * VM arguments in the Run Configuration.
 * 
 * @author Jiunn Haur Lim
 */
public class Sum {
  
  private static Core core;
  private static DirectedGraph<ScalarVertex, DefaultWeightedEdge> graph;

  public static void main(String[] args) throws Core.CoreException {
    
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    
    core = new Core();
    // JGraphT
    graph = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    int[] hello = {1, 2, 3, 4, 5};
    
    for (int i=0; i<hello.length; i++)
      graph.addVertex(new ScalarVertex(i, hello[i]));
    
    core.setGraph(graph);
    core.addAggregator("agg", new SumAggregator(0), 0);
    core.aggregateNow("agg");
    
    // don't forget to destroy!
    core.destroy();
    
  }
  
  private static class SumAggregator extends Aggregator<ScalarVertex, SumAggregator> {

    private int sum;
    
    public SumAggregator(int init){
      sum = init;
    }

    @Override
    protected void add(SumAggregator agg) {
      sum += agg.sum;
    }

    @Override
    protected SumAggregator clone() {
      return new SumAggregator(sum);
    }

    @Override
    protected void exec(Context context, ScalarVertex vertex) {
      sum += vertex.value();
    }

    @Override
    protected void finalize(Context context) {
      System.out.println("Sum: " + sum);
    }
    
  }
  
}