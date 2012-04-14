package org.graphlab.demo;

import java.io.IOException;
import java.util.Set;
import java.util.TreeSet;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.graphlab.util.GraphLoader;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Naive Greedy Graph Coloring Algorithm
 * 
 * <p>Demonstrates GraphLab Java. To run this class,
 * provide a path to a tsv file in the program arguments. Some examples are
 * available in the <tt>java_jni/test-graphs</tt> directory.</p>
 * 
 * <p>To run this demo, use</p>
<pre>
  ant Coloring -Dfile=&lt;path to tsv file&gt;
</pre>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class Coloring {

  private static final Logger logger = Logger.getLogger(Coloring.class);
  
  public static void main (String[] args){
    
    initLogger();
    
    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];

    // initialize graphlab core
    final Core core;
    try { core = new Core(); } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      return;
    }

    // construct graph
    final DirectedGraph<ScalarVertex, DefaultWeightedEdge> graph;
    try { graph = constructGraph(filename); } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      core.destroy();    // cleanup
      return;
    }
    
    for (ScalarVertex v : graph.vertexSet()) {
      v.setValue(-1);    // -1 means uncolored
    }
    
    // execute graph updates
    core.setGraph(graph);
    core.scheduleAll(new ColoringUpdater(graph));
    logger.info("Took " + core.start() + " seconds.");
    
    // print results and cleanup
    printResults (graph);
    core.destroy();
    
    return;

  }

  /**
   * Initializes the logger.
   */
  private static void initLogger() {
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
    logger.setLevel(Level.ALL);
  }
  
  /**
   * Prints results of algorithm.
   */
  private static void printResults (Graph<ScalarVertex, ?> graph){
    logger.info("----------------- Results -----------------");
    logger.info("ID : Color");
    for (ScalarVertex v : graph.vertexSet())
      logger.info(v.id() + " : " + (int) v.value());
  }
  
  /**
   * Checks that required input parameters are available and valid. Prints
   * instructions if not all parameters were valid.
   */
  private static boolean checkParams(String[] args) {
    if (args.length == 1) return true;
    System.out.println("Please provide filename.");
    System.out.println("Usage: ant Coloring -Dfile=path/to/tsv/file");
    return false;
  }

  private static DirectedGraph<ScalarVertex, DefaultWeightedEdge>
    constructGraph(String filename) throws IOException {

    DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    return graph;

  }
  
  private static class ColoringUpdater
    extends Updater<ScalarVertex, DefaultWeightedEdge, ColoringUpdater> {

    private DirectedGraph<ScalarVertex, DefaultWeightedEdge> mGraph;

    public ColoringUpdater(DirectedGraph<ScalarVertex, DefaultWeightedEdge> g) {
      this.mGraph = g;
    }

    @Override
    public void update(Context context, ScalarVertex vertex) {

      Set<Integer> neighborColors = new TreeSet<Integer>();
      int color;
      
      // collect neighbor colors
      for (DefaultWeightedEdge edge : mGraph.incomingEdgesOf(vertex)){
        // we will just use value as color
        color = (int) mGraph.getEdgeSource(edge).value();
        if (-1 == color) continue;
        neighborColors.add(color);
      }
      
      for (DefaultWeightedEdge edge : mGraph.outgoingEdgesOf(vertex)){
        color = (int) mGraph.getEdgeTarget(edge).value();
        if (-1 == color) continue;
        neighborColors.add(color);
      }
      
      // find a unique color
      for (color=0; ; color++){
        if (!neighborColors.contains(color)) break;
      }
      
      vertex.setValue(color);

    }

    @Override
    protected ColoringUpdater clone() {
      return new ColoringUpdater(mGraph);
    }

  }
  
}
