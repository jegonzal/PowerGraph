package org.graphlab.demo;

import java.io.IOException;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.logging.StreamHandler;

import org.graphlab.Core;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.graphlab.util.GraphLoader;

/**
 * ShortestPath algorithm. Demonstrates GraphLab over JNI.
 * @author Jiunn Haur Lim
 */
public class ShortestPath {

	public static void main(String[] args){
		
		if (args.length != 1){
			System.out.println ("Please provide filename.");
			return;
		}
		
		// setup logging
		initLogging();

		// init core
		final Core<SparseGraph<ScalarVertex, ScalarEdge>> c;
		try {
			c = new Core<SparseGraph<ScalarVertex, ScalarEdge>>();
		} catch (Exception e) {
			System.out.println ("Unable to initialize core.");
			return;
		}
		
		final SparseGraph<ScalarVertex, ScalarEdge> g = constructGraph(args[0]);
		if (null == g){
			System.out.println ("Unable to construct graph.");
			c.destroy();
			return;
		}
		
		Updater shortestPathUpdater = new Updater(){

			@Override
			public void update(int vertex_id) {
				
				ScalarVertex vertex = g.getVertex(vertex_id);
				
				for (ScalarEdge edge : g.incomingEdges(vertex.id())){
					ScalarVertex neighbor = g.getVertex(edge.source());
					vertex.setValue ((float) Math.min(vertex.value(), neighbor.value() + edge.weight()));
				}
				
			    // reschedule any affected neighbors
			    for (ScalarEdge edge : g.outgoingEdges(vertex_id)) {
			      ScalarVertex neighbor = g.getVertex(edge.target());
			      if (neighbor.value() > (vertex.value() + edge.weight())) 
			        c.schedule (edge.target(), this);    
			    }

			}
			
		};
		
		for (ScalarVertex v : g.vertices()){
			v.setValue(Integer.MAX_VALUE);
		}
		
		// start from root
		ScalarVertex root = g.getVertex(0);
		root.setValue(0);
		
		c.setGraph(g);
		c.schedule(root.id(), shortestPathUpdater);
		c.start();
		
		// destroy core
		c.destroy();
		
		System.out.println ("Shortest path from root to 7 was: " + g.getVertex(7).value());
		
	}
	
	private static void initLogging (){
		Handler handler = new StreamHandler(System.out, new SimpleFormatter());
		handler.setLevel(Level.ALL);
		Logger.getLogger(Core.TAG).addHandler(handler);
		Logger.getLogger(Core.TAG).setLevel(Level.ALL);
	}
	
	private static SparseGraph<ScalarVertex, ScalarEdge> constructGraph (String filename){
		
		SparseGraph<ScalarVertex, ScalarEdge> graph = new SparseGraph<ScalarVertex, ScalarEdge>();
		
		try {
			GraphLoader.loadGraphFromTsvFile(graph, filename);
		}catch (IOException e){
			System.out.println (e.getMessage());
			System.out.println ("Unable to load from file: " + filename);
			return null;
		}
		
		return graph;
		
	}
	
}
