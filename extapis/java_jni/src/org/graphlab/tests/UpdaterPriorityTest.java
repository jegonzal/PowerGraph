package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.util.concurrent.atomic.AtomicInteger;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Scheduler;
import org.graphlab.Updater;
import org.graphlab.data.Graph;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class UpdaterPriorityTest {

  private Core<Graph<ScalarVertex, ScalarEdge>> core;
  
  @Before
  public void setUp() throws Exception {
    // init logging
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    // create core
    core = new Core<Graph<ScalarVertex, ScalarEdge>>();
  }
  
  @Test
  public void testPriority(){
    
    final AtomicInteger flag = new AtomicInteger(0);
    
    // load graph
    final SparseGraph<ScalarVertex, ScalarEdge> graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(0));
    
    core.setSchedulerType(Scheduler.PRIORITY);
    core.setNCpus(1);
    core.setGraph(graph);
    
    core.schedule(0, new Updater(core){
      @Override
      public void update(Context context, int vertexId) {
        synchronized(flag){
          assertEquals(flag.getAndIncrement(), 1);
        }
      }
      @Override
      public double priority(){ return 6; }
    });
    
    core.schedule(1, new Updater(core){
      @Override
      public void update(Context context, int vertexId) {
        synchronized(flag){
          // this should be incremented first
          assertEquals(flag.getAndIncrement(), 0);
        }
      }
      @Override
      public double priority(){ return 10; }
    });
    
    core.start();
    
  }

  @After
  public void tearDown() throws Exception {
    core.destroy();
  }

}
