package org.graphlab.toolkits.matrix;

import java.io.FileReader;
import java.io.IOException;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import cern.colt.matrix.io.MatrixSize;
import cern.colt.matrix.io.MatrixVectorReader;

/**
 * Loads a matrix from a Matrix Market file and creates the corresponding graph.
 * @author Jiunn Haur Lim
 */
public class MatrixLoader {

  /**
   * Constructs graph from a MM file (assuming general coordinate format with real values.)
   * @param filename
   * @return undirected graph with vector vertices and scalar edges
   * @throws IOException
   */
  public static
    SimpleWeightedGraph<VectorVertex, DefaultWeightedEdge>
    loadGraph(String filename) throws IOException {
    
    MatrixVectorReader reader = new MatrixVectorReader(new FileReader(filename));
    MatrixSize size = reader.readMatrixSize(reader.readMatrixInfo());
    
    // create graph
    SimpleWeightedGraph<VectorVertex, DefaultWeightedEdge> graph
    = new SimpleWeightedGraph<VectorVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    
    // initialize graph size
    for (int i=0; i<size.numRows() + size.numColumns(); i++){
      graph.addVertex(new VectorVertex(i));
    }
    
    // iterate through file entries and construct graph
    int[] row = new int[1];
    int[] col = new int[1];
    double[] data = new double[1];
    for (int i=0; i<size.numEntries(); i++){
      
      reader.readCoordinate(row, col, data);
      
      // hopefully JGraphT uses their ID's to get the correct vertices
      VectorVertex source = new VectorVertex(row[0]);
      VectorVertex target = new VectorVertex(size.numRows()+col[0]);
      
      DefaultWeightedEdge edge = graph.addEdge(source, target);
      graph.setEdgeWeight(edge, data[0]);
      
    }
    
    return graph;
    
  }
  
}
