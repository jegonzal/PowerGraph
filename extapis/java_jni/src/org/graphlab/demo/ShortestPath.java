package org.graphlab.demo;

import java.io.IOException;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.Updater;
import org.graphlab.data.ScalarEdge;
import org.graphlab.data.ScalarVertex;
import org.graphlab.data.SparseGraph;
import org.graphlab.util.GraphLoader;

/**
 * ShortestPath algorithm.
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
 *  org.graphlab.demo.ShortestPath test-graphs/toy.tsv
 * </pre>
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ShortestPath {

  private static final Logger logger = Logger.getLogger(ShortestPath.class);

  public static void main(String[] args) {

    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
    logger.setLevel(Level.ALL);

    logger.trace("Main method in " + ShortestPath.class.getCanonicalName()
        + " started.");

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
      v.setValue(Integer.MAX_VALUE);
    }

    // start from root
    ScalarVertex root = g.getVertex(0);
    root.setValue(0);

    c.setGraph(g);
    c.schedule(root.id(), new ShortestPathUpdater(c, g));

    logger.trace("Running graphlab ...");
    c.start();

    // destroy core
    logger.trace("Destroying core ...");
    c.destroy();

    logger.info("Shortest path from root to 7 was: " + g.getVertex(7).value());

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
  
  private static class ShortestPathUpdater extends Updater {

    private SparseGraph<ScalarVertex, ScalarEdge> g;

    public ShortestPathUpdater(Core<?> core, SparseGraph<ScalarVertex, ScalarEdge> g) {
      super(core);
      this.g = g;
    }

    @Override
    public void update(Context context, int vertexId) {

      ScalarVertex vertex = g.getVertex(vertexId);

      // find shortest known distance into this node
      for (ScalarEdge edge : g.incomingEdges(vertex.id())) {
        ScalarVertex neighbor = g.getVertex(edge.source());
        vertex.setValue((float) Math.min(vertex.value(), neighbor.value()
            + edge.weight()));
      }

      // reschedule any affected neighbors
      for (ScalarEdge edge : g.outgoingEdges(vertexId)) {
        ScalarVertex neighbor = g.getVertex(edge.target());
        if (neighbor.value() > (vertex.value() + edge.weight())) {
          context.schedule(neighbor.id(), this);
        }
      }

    }

  }

}
