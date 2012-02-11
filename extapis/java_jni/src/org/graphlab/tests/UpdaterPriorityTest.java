package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Scheduler;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class UpdaterPriorityTest {

  private Core core;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    core = new Core();
  }
  
  @Test
  public void testPriority(){
    
    final AtomicInteger flag = new AtomicInteger(0);
    
    // load graph
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(1));
    
    core.setSchedulerType(Scheduler.PRIORITY);
    core.setNCpus(1);
    core.setGraph(graph);
    
    Iterator<ScalarVertex> it = graph.vertexSet().iterator();
    
    core.schedule(it.next(), new Updater<ScalarVertex>(){
      @Override
      public void update(Context context, ScalarVertex vertex) {
        synchronized(flag){
          assertEquals(flag.getAndIncrement(), 1);
        }
      }
      @Override
      public double priority(){ return 6; }
      @Override
      protected Updater<ScalarVertex> clone() {
        return this;
      }
    });
    
    core.schedule(it.next(), new Updater<ScalarVertex>(){
      @Override
      public void update(Context context, ScalarVertex v) {
        synchronized(flag){
          // this should be incremented first
          assertEquals(flag.getAndIncrement(), 0);
        }
      }
      @Override
      public double priority(){ return 10; }
      @Override
      protected Updater<ScalarVertex> clone() {
        return this;
      }
    });
    
    core.start();
    
  }

  @After
  public void tearDown() throws Exception {
    core.destroy();
  }

}
