package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.Context;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.Scope;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.junit.Before;
import org.junit.Test;

/**
 * Test cases that configure the core (via {@link Core#setNCpus(long)},
 * {@link Core#setSchedulerType(org.graphlab.Scheduler)}, and
 * {@link Core#setScopeType(Scope)}) after the core has already started.
 * 
 * These are weak test cases, which merely check that the engine doesn't crash
 * on a very simple dual vertex graph for each configuration. Manual checks are
 * recommended to ensure that the engine options (from the stdout) corresponds
 * to the test cases.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class BeforeStartConfigurationTest {

  private Core<SparseGraph<ScalarVertex, ScalarEdge>> core;
  private SparseGraph<ScalarVertex, ScalarEdge> graph;
  
  private static final Logger logger
    = Logger.getLogger(AfterStartConfigurationTest.class);

  @Before
  public void setUp() throws Exception {
    
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    
    // initialize graph to contain a single vertex of 0
    graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(0));
    
  }

  @Test
  public void testSetScheduler() throws CoreException{
    
    logger.debug("------ testSetScheduler BEGIN ------");
    
    for (Scheduler scheduler : Scheduler.values()){ 
      
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(scheduler);
      
      core = new Core<SparseGraph<ScalarVertex, ScalarEdge>>(config);
      core.setGraph(graph);
      
      // updater will set value to 1
      core.schedule(0, new Updater(core){
        @Override
        public void update(Context context, int vertexId) {
          graph.getVertex(vertexId).setValue(1);
        }
      });
      
      logger.debug("Expect scheduler: " + scheduler.toString());
      core.start();
      assertEquals (graph.getVertex(0).value(), 1.0, 0);
      
      core.destroy();
      
    }
    
    logger.debug("------ testSetScheduler BEGIN ------");
    
  }
  
  @Test
  public void testSetScopeType() throws CoreException{
    
    logger.debug("------ testSetScope BEGIN ------");
    
    for (Scope scope : Scope.values()){ 
      
      CoreConfiguration config = new CoreConfiguration();
      config.setScopeType(scope);
      
      core = new Core<SparseGraph<ScalarVertex, ScalarEdge>>(config);
      core.setGraph(graph);
      
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
      
      core.destroy();
      
    }
    
    logger.debug("------ testSetScope BEGIN ------");
    
  }
  
  @Test
  public void testSetNCpus() throws CoreException{
    
    logger.debug("------ testSetScope BEGIN ------");
    
    for (int n=1; n<=2; n++){ 
      
      CoreConfiguration config = new CoreConfiguration();
      config.setNCpus(n);
      
      core = new Core<SparseGraph<ScalarVertex, ScalarEdge>>(config);
      core.setGraph(graph);
      
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
      
      core.destroy();
      
    }
    
    logger.debug("------ testSetScope BEGIN ------");
    
  }

}
