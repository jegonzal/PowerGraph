package org.graphlab.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;

/**
 * Utility functions for loading graphs from various file formats.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class GraphLoader {

	/**
	 * Load a graph file specified in the format:
	 * 
	 * <code>
	 * source_id <tab> target_id
	 * source_id <tab> target_id
	 * source_id <tab> target_id
	 * </code>
	 * 
	 * The file should not contain repeated edges, and should not contain edge
	 * weights.
	 * 
	 * @param graph
	 *            the graph to mutate
	 * @param filename
	 *            the file to read from
	 */
	public static void loadGraphFromTsvFile(
			MutableGraph<ScalarVertex, ScalarEdge> graph, String filename)
			throws IOException {

		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String input;
		int lineNumber = 1;

		// TODO: update classics/shortest_path and then update this to include edge weights
		
		while (null != (input = reader.readLine())) {

			String[] tokens = input.split(" ");
			if (tokens.length != 2)
				throw new IOException("Line " + lineNumber + " did not have exactly two tokens.");

			int source = Integer.parseInt(tokens[0]);
			int target = Integer.parseInt(tokens[1]);

			// ensure that the number of vertices is correct
			if (source >= graph.size() || target >= graph.size())
				resizeGraph (graph, Math.max(source, target) + 1);
			if (source != target) {
				graph.addEdge(new ScalarEdge(source, target, 1));
			}

		}

		System.out.println("Finished loading graph with: " + graph.size()
				+ " vertices.");

	}
	
	private static void resizeGraph (MutableGraph<ScalarVertex, ScalarEdge> graph, int count){
		
		int difference = count - graph.size();
		while (difference > 0){
			graph.addVertex(new ScalarVertex(0));
			difference--;
		}
		
	}

}
