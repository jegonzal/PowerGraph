package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.Graph;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class AddUpdaterTest {

  private Core<Graph<ScalarVertex, ScalarEdge>> core;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    core = new Core<Graph<ScalarVertex, ScalarEdge>>();
  }

  /**
   * Schedules the same vertex twice and checks that the updater counter is incremented
   * correctly.
   */
  @Test
  public void testSimpleAddUpdater(){
    
    // load graph
    final SparseGraph<ScalarVertex, ScalarEdge> graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(1));

    core.setGraph(graph);
    core.schedule(0, new CounterUpdater(core, graph));
    core.start();
    
    // check shortest paths
    assertEquals (null, 4, graph.getVertex(1).value(), 0);
    
  }
  
  @After
  public void tearDown() throws Exception {
    core.destroy();
  }
  
  private static class CounterUpdater extends Updater {

    private Graph<ScalarVertex, ScalarEdge> mGraph;
    
    public CounterUpdater(Core<Graph<ScalarVertex, ScalarEdge>> core,
                          Graph<ScalarVertex, ScalarEdge> graph) {
      super(core);
      mGraph = graph;
    }

    private int counter = 1;

    @Override
    public void update(Context context, int vertexId) {
      // record counter value in vertex
      mGraph.getVertex(vertexId).setValue(counter);
      if (0 == vertexId){
        // schedule one three times
        context.schedule(1, this);      // count = 1
        context.schedule(1, this);      // count = 1 + 1 = 2
        context.schedule(1, this);      // count = 2 + 2 = 4
      }
    }

    @Override
    public void add(Updater updater) {
      if (!(updater instanceof CounterUpdater)) return;
      CounterUpdater cu = (CounterUpdater) updater;
      counter += cu.counter;
    }

  }

}
