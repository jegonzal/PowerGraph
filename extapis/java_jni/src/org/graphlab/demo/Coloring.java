package org.graphlab.demo;

import java.io.IOException;
import java.util.Set;
import java.util.TreeSet;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Scope;
import org.graphlab.Core.CoreException;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.graphlab.util.GraphLoader;

/**
 * Naive Greedy Graph Coloring Algorithm
 * 
 * <p>Demonstrates GraphLab over JNI. To run this class,
 * provide a path to a tsv file in the program arguments. Some examples are
 * available in the <tt>java_jni/test-graphs</tt> directory.</p>
 * 
 * <p>For example:
 * <pre>
 * cd extapis/java_jni
 * java \
 *  -classpath bin:lib/log4j-1.2.16.jar \
 *  -Djava.library.path=../../release/src/graphlab/jni \
 *  org.graphlab.demo.Coloring test-graphs/toy.tsv
 * </pre>
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class Coloring {

  private static final Logger logger = Logger.getLogger(Coloring.class);
  
  public static void main (String[] args){
    
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
    logger.setLevel(Level.ALL);
    
    // check arguments
    if (!checkParams(args)) {
      logger.trace("Exiting main method.");
      return;
    }

    String filename = args[0];
    logger.info("Graph file: " + filename);

    // initialize graphlab core
    final Core<SparseGraph<ScalarVertex, ScalarEdge>> c;
    try {
      logger.trace("Initializing GraphLab core ...");
      c = new Core<SparseGraph<ScalarVertex, ScalarEdge>>();
    } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      logger.trace("Exiting main method.");
      return;
    }

    // construct graph
    logger.trace("Constructing graph from " + filename + " ...");
    final SparseGraph<ScalarVertex, ScalarEdge> g;
    try {
      g = constructGraph(filename);
    } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      c.destroy();
      logger.trace("Exiting main method.");
      return;
    }
    
    for (ScalarVertex v : g.vertices()) {
      v.setValue(-1);    // -1 means uncolored
    }
    
    // execute graph updates
    c.setScopeType(Scope.EDGE);
    c.setGraph(g);
    c.scheduleAll(new ColoringUpdater(g));
    logger.trace("Running graphlab ...");
    logger.info("Took " + c.start() + " seconds.");
    
    // print results
    printResults (g);

  }
  
  private static void printResults (SparseGraph<ScalarVertex, ScalarEdge> g){
    
    logger.info("----------------- Results -----------------");
    logger.info("ID : Color");
    for (ScalarVertex v : g.vertices()){
      logger.info(v.id() + " : " + (int) v.value());
    }
    
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
          + ShortestPath.class.getCanonicalName() + " path/to/tsv/file");
      return false;
    }

    return true;

  }

  private static SparseGraph<ScalarVertex, ScalarEdge> constructGraph(
      String filename) throws IOException {

    SparseGraph<ScalarVertex, ScalarEdge> graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    return graph;

  }
  
  private static class ColoringUpdater extends Updater {

    private SparseGraph<ScalarVertex, ScalarEdge> g;

    public ColoringUpdater(SparseGraph<ScalarVertex, ScalarEdge> g) {
      this.g = g;
    }

    @Override
    public void update(Context context, int vertexId) {

      ScalarVertex vertex = g.getVertex(vertexId);
      Set<Integer> neighborColors = new TreeSet<Integer>();
      int color;
      
      // collect neighbor colors
      for (ScalarEdge edge : g.incomingEdges(vertexId)){
        // we will just use value as color
        color = (int) g.getVertex(edge.source()).value();
        if (-1 == color) continue;
        neighborColors.add(color);
      }
      
      for (ScalarEdge edge : g.outgoingEdges(vertexId)){
        color = (int) g.getVertex(edge.target()).value();
        if (-1 == color) continue;
        neighborColors.add(color);
      }
      
      // find a unique color
      for (color=0; ; color++){
        if (!neighborColors.contains(color)) break;
      }
      
      vertex.setValue(color);

    }

  }
  
}
