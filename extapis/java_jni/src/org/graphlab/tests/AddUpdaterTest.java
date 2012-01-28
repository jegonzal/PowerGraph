package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.util.Iterator;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class AddUpdaterTest {

  private Core core;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    core = new Core();
  }

  /**
   * Schedules the same vertex twice and checks that the updater counter is incremented
   * correctly.
   */
  @Test
  public void testSimpleAddUpdater(){
    
    // load graph
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(1));

    Iterator<ScalarVertex> it = graph.vertexSet().iterator();
    
    core.setGraph(graph);
    core.schedule(it.next(), new CounterUpdater(graph));
    core.start();
    
    // check counter value
    assertEquals (null, 4, it.next().value(), 0);
    
  }
  
  @After
  public void tearDown() throws Exception {
    core.destroy();
  }
  
  private static class CounterUpdater extends Updater<ScalarVertex> {

    private Graph<ScalarVertex, DefaultWeightedEdge> mGraph;
    
    public CounterUpdater(Graph<ScalarVertex, DefaultWeightedEdge> graph) {
      mGraph = graph;
    }

    private int counter = 1;

    @Override
    public void update(Context context, ScalarVertex vertex) {
      
      Iterator<ScalarVertex> it = mGraph.vertexSet().iterator();
      it.next();
      ScalarVertex one = it.next();
      
      // record counter value in vertex
      vertex.setValue(counter);
      if (0 == vertex.id()){
        // schedule one three times
        context.schedule(one, this);      // count = 1
        context.schedule(one, this);      // count = 1 + 1 = 2
        context.schedule(one, this);      // count = 2 + 2 = 4
      }
    }

    @Override
    public void add(Updater<ScalarVertex> updater) {
      if (!(updater instanceof CounterUpdater)) return;
      CounterUpdater cu = (CounterUpdater) updater;
      counter += cu.counter;
    }

  }

}
