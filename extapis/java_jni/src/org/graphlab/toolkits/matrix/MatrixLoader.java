package org.graphlab.toolkits.matrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import org.graphlab.data.Vertex;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Loads a matrix from a MatrixMarket file and constructs the corresponding graph.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class MatrixLoader {

  /**
   * Constructs graph from a TSV file.
   * 
   * Format: user <tab> item <tab> rating
   * 
   * @param graph       the graph to construct
   * @param numRows     number of rows
   * @param numCols     number of columns
   * @param filename    name of file containing the matrix
   * @throws IOException
   *            if file could not be read
   * @throws IllegalAccessException 
   *            if vertices could not be instantiated
   * @throws InstantiationException 
   *            if vertices could not be instantiated
   */
  public static <V extends Vertex> void
    loadGraphFromTsv(
        WeightedGraph<V, DefaultWeightedEdge> graph,
        Class<V> vertexClass,
        int numRows,
        int numCols,
        String filename)
  throws IOException, InstantiationException, IllegalAccessException{
    
    if (null == graph || null == vertexClass || null == filename)
      throw new NullPointerException("graph, vertexClass, and filename cannot be null.");
    
    // read from file
    BufferedReader reader = new BufferedReader(new FileReader(filename));
    Map<Integer, V> vertices = new HashMap<Integer, V>();
    
    String input;
    int lineNumber = 1;
    
    // read line by line
    while (null != (input = reader.readLine())) {

      // break up into source - target - weight
      String[] tokens = input.trim().split("\\s+");
      if (tokens.length < 3)
        throw new IOException("Error parsing line " + lineNumber);

      int source = Integer.parseInt(tokens[0]);
      int target = Integer.parseInt(tokens[1]) + numRows;
      double weight = Double.parseDouble(tokens[2]);
      
      // create vertices if encountering them for the first time
      V srcV = vertices.get(source);
      if (null == srcV){
        srcV = vertexClass.newInstance();
        srcV.setId(source);
        graph.addVertex(srcV);
        vertices.put(source, srcV);
      }
      
      // create target vertex
      V trgtV = vertices.get(target);
      if (null == trgtV){
          trgtV = vertexClass.newInstance();
          trgtV.setId(target);
          graph.addVertex(trgtV);
          vertices.put(target, trgtV);
      }
      
      DefaultWeightedEdge edge = graph.addEdge(srcV, trgtV);
      if (null != edge)
        graph.setEdgeWeight(edge, weight);

    }
    
    // add the remaining vertices
    for (int i=0; i<numRows + numCols; i++){
      
      if (vertices.containsKey(i)) continue;
      
      V vertex = vertexClass.newInstance();
      vertex.setId(i);
      graph.addVertex(vertex);
      
    }
    
  }
  
  /**
   * Constructs graph from an MM file (assuming general coordinate format with real values.)
   * 
   * @param graph       the graph to construct
   * @param filename    name of file containing the matrix
   * @throws IOException
   *            if file could not be read
   * @throws IllegalAccessException 
   *            if vertices could not be instantiated
   * @throws InstantiationException 
   *            if vertices could not be instantiated
   */
  public static <V extends Vertex> void
    loadGraphFromMM(
        WeightedGraph<V, DefaultWeightedEdge> graph,
        Class<V> vertexClass,
        String filename)
  throws IOException, InstantiationException, IllegalAccessException {
    
    if (null == graph || null == vertexClass || null == filename)
      throw new NullPointerException("graph, vertexClass, and filename cannot be null.");
    
    // read matrix metadata
    MatrixVectorReader reader = new MatrixVectorReader(new FileReader(filename));
    MatrixSize size = reader.readMatrixSize(reader.readMatrixInfo());
    
    // I really need this - hashcode trick doesn't work
    Map<Integer, V> vertices = new HashMap<Integer, V>();
    
    // iterate through file entries and construct graph
    int[] row = new int[1];
    int[] col = new int[1];
    double[] data = new double[1];
    for (int i=0; i<size.numEntries(); i++){
      
      reader.readCoordinate(row, col, data);
      
      // add vertex if it does not already exist
      V source = vertices.get(row[0]);
      if (null == source){
        source = vertexClass.newInstance();
        source.setId(row[0]);
        graph.addVertex(source);
        vertices.put(source.id(), source);
      }
      
      // add vertex if it does not already exist
      int targetID = size.numRows() + col[0];
      V target = vertices.get(targetID);
      if (null == target){
         target = vertexClass.newInstance();
         target.setId(targetID);
         graph.addVertex(target);
         vertices.put(targetID, target);
      }
      
      DefaultWeightedEdge edge = graph.addEdge(source, target);
      graph.setEdgeWeight(edge, data[0]);
      
    }
    
    // add the remaining vertices
    for (int i=0; i<size.numRows() + size.numColumns(); i++){
      
      if (vertices.containsKey(i)) continue;
      
      V vertex = vertexClass.newInstance();
      vertex.setId(i);
      graph.addVertex(vertex);
      
    }
    
  }
  
}
