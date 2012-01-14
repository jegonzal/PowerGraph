package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Scheduler;
import org.graphlab.Core.CoreException;
import org.graphlab.Scope;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Test cases that reconfigure the core (via {@link Core#setNCpus(long)}, 
 * {@link Core#setSchedulerType(org.graphlab.Scheduler)}, and
 * {@link Core#setScopeType(Scope)}) after the core has already started.
 * 
 * These are weak test cases, which merely check that the engine doesn't
 * crash on a very simple single vertex graph for each configuration. Manual
 * checks are recommended to ensure that the engine options (from the stdout)
 * corresponds to the test cases.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class AfterStartConfigurationTest {

  private static final Logger logger
    = Logger.getLogger(AfterStartConfigurationTest.class);

  private Core<SparseGraph<ScalarVertex, ScalarEdge>> core;
  private SparseGraph<ScalarVertex, ScalarEdge> graph;
  
  @Before
  public void setUp() throws CoreException {

    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    
    // initialize core
    core = new Core<SparseGraph<ScalarVertex, ScalarEdge>>();
    
    // initialize graph to contain a single vertex of 0
    graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(0));
    
    core.setGraph(graph);
    core.start();
    
  }
  
  @Test
  public void testSetScopeAfterStart(){
    
    logger.debug("------ testSetScopeAfterStart BEGIN ------");
    
    for (Scope scope : Scope.values()){
      
      // initialize vertex to 0
      graph.getVertex(0).setValue(0);
      
      core.setScopeType(scope);
      
      // updater will set value to 1
      core.schedule(0, new Updater(core){
        @Override
        public void update(Context context, int vertexId) {
          graph.getVertex(vertexId).setValue(1);
        }
      });
      
      logger.debug("Expect scope: " + scope.toString());
      core.start();
      
      assertEquals (graph.getVertex(0).value(), 1.0, 0);
      
    }
    
    logger.debug("------ testSetScopeAfterStart END --------");
    
  }
  
  @Test
  public void testSetSchedulerAfterStart(){
    
    logger.debug("------ testSetSchedulerAfterStart BEGIN ------");
    
    for (Scheduler scheduler : Scheduler.values()){
      
      // initialize vertex to 0
      graph.getVertex(0).setValue(0);
      
      core.setSchedulerType(scheduler);
      
      // updater will set value to 1
      core.schedule(0, new Updater(core){
        @Override
        public void update(Context context, int vertexId) {
          graph.getVertex(vertexId).setValue(1);
        }
      });
      
      logger.debug("Expect scheduler: " + scheduler.type());
      core.start();
      
      assertEquals (graph.getVertex(0).value(), 1.0, 0);
      
    }
    
    logger.debug("------ testSetSchedulerAfterStart END --------");
    
  }
  
  @Test
  public void testSetNCpusAfterStart(){
   
    logger.debug("------ testSetNCpusAfterStart BEGIN ------");
    
    for (int n=1; n<=2; n++){
      
      // initialize vertex to 0
      graph.getVertex(0).setValue(0);
      
      core.setNCpus(n);
      
      // updater will set value to 1
      core.schedule(0, new Updater(core){
        @Override
        public void update(Context context, int vertexId) {
          graph.getVertex(vertexId).setValue(1);
        }
      });
      
      logger.debug("Expect ncpus: " + n);
      core.start();
      
      assertEquals (graph.getVertex(0).value(), 1.0, 0);
      
    }
    
    logger.debug("------ testSetNCpusAfterStart END --------");
    
  }

  @After
  public void tearDown() throws Exception {
    core.destroy();
  }

}
