package org.graphlab.toolkits.matrix;

import java.io.FileReader;
import java.io.IOException;

import org.graphlab.data.Vertex;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

import cern.colt.matrix.io.MatrixSize;
import cern.colt.matrix.io.MatrixVectorReader;

/**
 * Loads a matrix from a MatrixMarket file and populates the corresponding graph.
 * @author Jiunn Haur Lim
 */
public class MatrixLoader {

  /**
   * Constructs graph from an MM file (assuming general coordinate format with real values.)
   * @param graph       the graph to construct
   * @param filename    name of file containing the matrix
   * @throws IOException
   * @throws IllegalAccessException 
   * @throws InstantiationException 
   */
  public static <V extends Vertex> void
    loadGraph(
        WeightedGraph<V, DefaultWeightedEdge> graph,
        Class<V> vertexClass,
        String filename)
  throws IOException, InstantiationException, IllegalAccessException {
    
    if (null == graph || null == vertexClass || null == filename)
      throw new NullPointerException("graph, vertexClass, and filename cannot be null.");
    
    // read matrix metadata
    MatrixVectorReader reader = new MatrixVectorReader(new FileReader(filename));
    MatrixSize size = reader.readMatrixSize(reader.readMatrixInfo());
    
    // initialize graph size
    for (int i=0; i<size.numRows() + size.numColumns(); i++){
      V vertex = vertexClass.newInstance();
      vertex.setId(i);
      graph.addVertex(vertex);
    }
    
    // iterate through file entries and construct graph
    int[] row = new int[1];
    int[] col = new int[1];
    double[] data = new double[1];
    for (int i=0; i<size.numEntries(); i++){
      
      reader.readCoordinate(row, col, data);
      
      // hopefully JGraphT uses their ID's to get the correct vertices
      V source = vertexClass.newInstance();
      source.setId(row[0]);
      V target = vertexClass.newInstance();
      target.setId(size.numRows()+col[0]);
      
      DefaultWeightedEdge edge = graph.addEdge(source, target);
      graph.setEdgeWeight(edge, data[0]);
      
    }
    
  }
  
}
