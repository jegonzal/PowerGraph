package org.graphlab.tests;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.graphlab.util.GraphLoader;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Simple deterministic test case for a shortest path algorithm using GraphLab.
 * 
 * @author Jiunn Haur Lim
 */
public class ShortestPathTest {

	private Core<SparseGraph<ScalarVertex, ScalarEdge>> c;
	
	@Before
	public void setUp() throws Exception {
		// init logging
		BasicConfigurator.configure();
		Logger.getLogger(Core.class).setLevel(Level.OFF);
		// create core
		c = new Core<SparseGraph<ScalarVertex, ScalarEdge>>();
	}

	@Test
	public void testSingleEdge() throws IOException {
		
		// load graph
		SparseGraph<ScalarVertex, ScalarEdge> g = new SparseGraph<ScalarVertex, ScalarEdge>();
		GraphLoader.loadGraphFromTsvFile(g, "test-graphs/one.tsv");
		
		// create updater
		Updater shortestPathUpdater = new ShortestPathUpdater(g);
		
		// initialize to infinity
		for (ScalarVertex v : g.vertices()) {
			v.setValue(Integer.MAX_VALUE);
		}

		// start from root
		ScalarVertex root = g.getVertex(0);
		root.setValue(0);

		c.setGraph(g);
		c.schedule(root.id(), shortestPathUpdater);
		c.start();
		
		// check shortest paths
		System.out.println (g);
		assertEquals (null, 10, g.getVertex(1).value(), 0);
		
	}
	
	// @Test
	public void testToyGraph() throws IOException {

		// load graph
		SparseGraph<ScalarVertex, ScalarEdge> g = new SparseGraph<ScalarVertex, ScalarEdge>();
		GraphLoader.loadGraphFromTsvFile(g, "test-graphs/toy.tsv");
		
		// create updater
		Updater shortestPathUpdater = new ShortestPathUpdater(g);
		
		// initialize to infinity
		for (ScalarVertex v : g.vertices()) {
			v.setValue(Integer.MAX_VALUE);
		}

		// start from root
		ScalarVertex root = g.getVertex(0);
		root.setValue(0);

		c.setGraph(g);
		c.schedule(root.id(), shortestPathUpdater);
		c.start();
		
		// check shortest paths
		assertEquals (null, 10, g.getVertex(2).value(), 0);
		assertEquals (null, 5, g.getVertex(3).value(), 0);
		assertEquals (null, 15, g.getVertex(4).value(), 0);
		assertEquals (null, 15, g.getVertex(5).value(), 0);
		assertEquals (null, 12, g.getVertex(6).value(), 0);
		assertEquals (null, 14, g.getVertex(7).value(), 0);

	}

	@After
	public void tearDown() throws Exception {
		c.destroy();
	}

	private class ShortestPathUpdater extends Updater {
	  
		private SparseGraph<ScalarVertex, ScalarEdge> g;

		public ShortestPathUpdater (SparseGraph<ScalarVertex, ScalarEdge> g) {
			this.g = g;
		}

		@Override
		public void update (long corePtr, long contextPtr, int vertexId){

			ScalarVertex vertex = g.getVertex(vertexId);
			System.out.println ("Updater called on vertex " + vertexId);

			// relax edges coming into this node
			for (ScalarEdge edge : g.incomingEdges(vertex.id())) {
				ScalarVertex neighbor = g.getVertex(edge.source());
				vertex.setValue
					((float) Math.min
						(vertex.value(),
						neighbor.value() + edge.weight()));
			}

			// reschedule any affected neighbors
			for (ScalarEdge edge : g.outgoingEdges(vertexId)) {
				ScalarVertex neighbor = g.getVertex(edge.target());
				System.out.println ("Checking vertex " + neighbor.id());
				if (neighbor.value() > (vertex.value() + edge.weight())){
					System.out.println ("rescheduling vertex " + neighbor.id());
					Context.getInstance().schedule(corePtr, contextPtr, neighbor.id(), this);
				}
			}

		}

	}

}
