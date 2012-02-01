package org.graphlab.tests;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class CoreTest {

  private Core mCore;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    mCore = new Core();
  }
  
  @Test
  public void testLastUpdateCount(){
    
    // create graph with 10 vertices
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    for (int i=0; i<10; i++) graph.addVertex(new ScalarVertex(i));
    
    // schedule a simple updater on all 10 vertices
    mCore.setGraph(graph);
    mCore.scheduleAll(new Updater<ScalarVertex>(){
      @Override
      public void update(Context context, ScalarVertex vertex){}
    });
    mCore.start();
    
    // check count
    assertEquals("Checking lastUpdateCount", 10, mCore.lastUpdateCount());
    
  }

  @After
  public void tearDown() throws Exception {
    mCore.destroy();
  }

}
