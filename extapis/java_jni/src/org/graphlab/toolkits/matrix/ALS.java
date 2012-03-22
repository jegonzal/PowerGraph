package org.graphlab.toolkits.matrix;

import java.io.IOException;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;


/**
 * Matrix factorization w. alternating least squares.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ALS {
  
  /**
   * Executes ALS matrix factorization on input file.
   * @param args
   *          command line arguments
   * @throws IOException 
   */
  public static void main(String[] args) throws IOException{
    
    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];
    
    SimpleWeightedGraph<VectorVertex, DefaultWeightedEdge>
      graph = MatrixLoader.loadGraph(filename);
    
  }
  
  /**
   * Checks that required input parameters are available and valid. Prints
   * instructions if not all parameters were valid.
   * 
   * @param args
   *          array of program arguments
   * @return true if parameters are OK; false otherwise
   */
  private static boolean checkParams(String[] args) {

    if (args.length != 1) {
      System.out.println("Please provide filename.");
      System.out.println("Usage: java -Djava.library.path=... "
          + ALS.class.getCanonicalName() + " path/to/matrix/file");
      return false;
    }

    return true;

  }

}
