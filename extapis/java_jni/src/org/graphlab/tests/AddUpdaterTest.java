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
    
    mCore.setGraph(graph);
    
    // schedule updater on vertex 0
    mCore.schedule(it.next(), new CounterUpdater(graph));
    mCore.start();
    
    // check counter value on one
    assertEquals ("Checking last counter value of updater.", 4, it.next().value(), 0);
    
  }
  
  @After
  public void tearDown() throws Exception {
    mCore.destroy();
  }
  
  /**
   * Counter Updater. When fused, counter values are added. When it
   * updates a vertex, it sets the value of the vertex with the current
   * value of its counter.
   * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
   */
  private static class CounterUpdater extends Updater<ScalarVertex> {

    private Graph<ScalarVertex, DefaultWeightedEdge> mGraph;
    
    public CounterUpdater(Graph<ScalarVertex, DefaultWeightedEdge> graph) {
      mGraph = graph;
    }

    private int counter = 1;

    @Override
    public void update(Context context, ScalarVertex vertex) {
      
      Iterator<ScalarVertex> it = mGraph.vertexSet().iterator();
      // skip vertex 0
      it.next();
      // get vertex 1
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

    /*
     * Fuse two updaters by adding count values
     * (non-Javadoc)
     * @see org.graphlab.Updater#add(org.graphlab.Updater)
     */
    @Override
    public void add(Updater<ScalarVertex> updater) {
      if (!(updater instanceof CounterUpdater)) return;
      CounterUpdater cu = (CounterUpdater) updater;
      counter += cu.counter;
    }

  }

}
