package org.graphlab.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;
import org.graphlab.data.ScalarVertex;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Utility functions for loading graphs from various file formats.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class GraphLoader {

	private static final Logger logger = Logger.getLogger(GraphLoader.class);
	
	/**
	 * Load a graph file specified in the format:
	 * 
	 * <code>
	 * source_id <tab> target_id <tab> weight
	 * source_id <tab> target_id <tab> weight
	 * source_id <tab> target_id <tab> weight
	 * </code>
	 * 
	 * The file should not contain repeated edges, and should not contain 
	 * self edges.
	 * 
	 * @param graph
	 *            the graph to construct
	 * @param filename
	 *            the file to read from
	 */
	public static void loadGraphFromTsvFile(
	    WeightedGraph<ScalarVertex, DefaultWeightedEdge> graph,
	    String filename)
			throws IOException {

		// read from file
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		// map from number to vertex
		Map<Integer, ScalarVertex> vertices = new HashMap<Integer, ScalarVertex>();
		
		String input;
		int lineNumber = 1;
		boolean warned = false;
		
		// read line by line
		while (null != (input = reader.readLine())) {

		  // break up into source - target - weight
			String[] tokens = input.trim().split("\\s+");
			if (tokens.length != 2 && tokens.length != 3)
				throw new IOException("Line " + lineNumber + " did not have exactly two/three tokens.");

			int source = Integer.parseInt(tokens[0]);
			int target = Integer.parseInt(tokens[1]);
			float weight = 1;
			if (3 == tokens.length) weight = Float.parseFloat(tokens[2]);

			// create vertices if encountering them for the first time
			ScalarVertex srcV = vertices.get(source);
			if (null == srcV){
			  srcV = new ScalarVertex(source);
			  vertices.put(source, srcV);
			  graph.addVertex(srcV);
			}
			
			ScalarVertex trgtV = vertices.get(target);
			if (null == trgtV){
			  trgtV = new ScalarVertex(target);
			  vertices.put(target, trgtV);
			  graph.addVertex(trgtV);
			}
			
			if (source != target) {
			  DefaultWeightedEdge edge = graph.addEdge(srcV, trgtV);
			  graph.setEdgeWeight(edge, weight);
			}else {
			  if (!warned){
			    logger.warn("Dropped self-edge for vertex " + source +
			                ". Subsequent warnings will be suppressed.");
			    warned = true;
			  }
			}

		}

		logger.trace ("Finished loading graph with: " + graph.vertexSet().size() + " vertices.");

	}

}
