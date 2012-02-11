package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Iterator;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.graphlab.util.GraphLoader;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Simple deterministic test case for a shortest path algorithm using GraphLab.
 * 
 * @author Jiunn Haur Lim
 */
public class ShortestPathTest {

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
	public void testSingleEdge() throws IOException {
		
		// load graph
		DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
		  = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		GraphLoader.loadGraphFromTsvFile(graph, "test-graphs/one.tsv");
		
		// initialize all distances to infinity
		for (ScalarVertex v : graph.vertexSet()) {
			v.setValue(Integer.MAX_VALUE);
		}

		// start from root, where distance=0
		ScalarVertex root = graph.vertexSet().iterator().next();
		root.setValue(0);

		core.setGraph(graph);
		Updater<ScalarVertex> shortestPathUpdater = new ShortestPathUpdater(graph);
		core.schedule(root, shortestPathUpdater);
		core.start();
		
		// check shortest paths
		Iterator<ScalarVertex> it = graph.vertexSet().iterator();
		it.next();  // get vertex 1
		assertEquals (null, 10, it.next().value(), 0);
		
	}
	
	@Test
	public void testToyGraph() throws IOException {

		// load graph
	  DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
	    = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		GraphLoader.loadGraphFromTsvFile(graph, "test-graphs/toy.tsv");
		
		// initialize all distances to infinity
		for (ScalarVertex v : graph.vertexSet()) {
			v.setValue(Integer.MAX_VALUE);
		}

		// start from root, where distance=0
		ScalarVertex root = graph.vertexSet().iterator().next();
		root.setValue(0);

		core.setGraph(graph);
		Updater<ScalarVertex> shortestPathUpdater = new ShortestPathUpdater(graph);
		core.schedule(root, shortestPathUpdater);
		core.start();
		
		// check shortest paths
		int[] expected = {0, 100, 10, 5, 15, 15, 12, 14};
		Iterator<ScalarVertex> it = graph.vertexSet().iterator();
		ScalarVertex vertex;
		for (int i=0; i<expected.length; i++){
		  vertex = it.next();
		  assertEquals(null, expected[i], vertex.value(), 0);
		}

	}

	@After
	public void tearDown() throws Exception {
		core.destroy();
	}

	private class ShortestPathUpdater extends Updater<ScalarVertex> {
	  
		private DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph;

		public ShortestPathUpdater
		  (DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> g) {
			this.graph = g;
		}

		@Override
		public void update (Context context, ScalarVertex vertex){

			// relax edges coming into this node
			for (DefaultWeightedEdge edge : graph.incomingEdgesOf(vertex)) {
				ScalarVertex neighbor = graph.getEdgeSource(edge);
				vertex.setValue
					((float)
					  Math.min(
					    vertex.value(),
						  neighbor.value() + graph.getEdgeWeight(edge)
					   )
					);
			}

			// reschedule any affected neighbors
			for (DefaultWeightedEdge edge : graph.outgoingEdgesOf(vertex)) {
				ScalarVertex neighbor = graph.getEdgeTarget(edge);
				if (neighbor.value() > (vertex.value() + graph.getEdgeWeight(edge))){
					context.schedule(neighbor, this);
				}
			}

		}

    @Override
    protected Updater<ScalarVertex> clone() {
      return new ShortestPathUpdater(graph);
    }

	}

}
