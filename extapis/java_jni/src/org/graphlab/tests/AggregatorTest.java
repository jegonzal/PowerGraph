package org.graphlab.tests;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Aggregator;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class AggregatorTest {

  private Core mCore;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    mCore = new Core();
  }
  
  /**
   * Simple Sum
   */
  @Test
  public void testAggregator(){
    
    // create graph with 10 vertices, each with values 0 to 9
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    for (int i=0; i<10; i++) graph.addVertex(new ScalarVertex(i, i));
    
    // run a simple aggregator on all 10 vertices
    mCore.setGraph(graph);
    mCore.addAggregator("aggregator", new SumAggregator(45), 0);
    mCore.aggregateNow("aggregator");
    
  }

  @After
  public void tearDown() throws Exception {
    mCore.destroy();
  }
  
  private class SumAggregator extends Aggregator<ScalarVertex>{

    private int mSum = 0;
    private final int mExpected;
    
    public SumAggregator(int expected){
      mExpected = expected;
    }
    
    
    @Override
    protected void exec(Context context, ScalarVertex vertex) {
      mSum += vertex.value();
    }

    @Override
    protected void add(Aggregator<ScalarVertex> aggregator) {
      if (aggregator instanceof SumAggregator)
         mSum += ((SumAggregator) aggregator).mSum; 
    }

    @Override
    protected void finalize(Context context) {
      // check sum
      assertEquals("Comparing sum against expected sum", mExpected, mSum);
    }
    
    @Override
    public Aggregator<ScalarVertex> clone(){
      return new SumAggregator(mExpected);
    }
    
  }

}
