package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.util.Iterator;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.Scope;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
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

  private Core core;
  private DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph;
  
  private static final Logger logger
    = Logger.getLogger(AfterStartConfigurationTest.class);

  @Before
  public void setUp() throws Exception {
    
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.OFF);
    
    // initialize graph to contain a single vertex of 0
    graph = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    graph.addVertex(new ScalarVertex(0));
    graph.addVertex(new ScalarVertex(0));
    
  }

  @Test
  public void testSetScheduler() throws CoreException{
    
    logger.debug("------ testSetScheduler BEGIN ------");
    
    for (Scheduler scheduler : Scheduler.values()){ 
      
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(scheduler);
      
      core = new Core(config);
      core.setGraph(graph);
      
      Iterator<ScalarVertex> it = graph.vertexSet().iterator();
      
      // updater will set value to 1
      core.schedule(it.next(), new Updater<ScalarVertex>(){
        @Override
        public void update(Context context, ScalarVertex vertex) {
          vertex.setValue(1);
        }
      });
      
      logger.debug("Expect scheduler: " + scheduler.toString());
      core.start();
      
      it = graph.vertexSet().iterator();
      assertEquals (it.next().value(), 1.0, 0);
      
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
      
      core = new Core(config);
      core.setGraph(graph);
      
      // updater will set value to 1
      Iterator<ScalarVertex> it = graph.vertexSet().iterator();
      core.schedule(it.next(), new Updater<ScalarVertex>(){
        @Override
        public void update(Context context, ScalarVertex vertex) {
          vertex.setValue(1);
        }
      });
      
      logger.debug("Expect scope: " + scope.toString());
      core.start();
      
      it = graph.vertexSet().iterator();
      assertEquals (it.next().value(), 1.0, 0);
      
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
      
      core = new Core(config);
      core.setGraph(graph);
      
      // updater will set value to 1
      Iterator<ScalarVertex> it = graph.vertexSet().iterator();
      core.schedule(it.next(), new Updater<ScalarVertex>(){
        @Override
        public void update(Context context, ScalarVertex vertex) {
          vertex.setValue(1);
        }
      });
      
      logger.debug("Expect ncpus: " + n);
      core.start();
      
      it = graph.vertexSet().iterator();
      assertEquals (it.next().value(), 1.0, 0);
      
      core.destroy();
      
    }
    
    logger.debug("------ testSetScope BEGIN ------");
    
  }

}
