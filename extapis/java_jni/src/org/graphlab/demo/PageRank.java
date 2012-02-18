package org.graphlab.demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.Updater;
import org.graphlab.data.ScalarVertex;
import org.graphlab.util.GraphLoader;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * PageRank Algorithm
 * 
 * <p>Demonstrates GraphLab over JNI.</p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class PageRank {

  private static final Logger logger = Logger.getLogger(PageRank.class);
  
  public static void main (String[] args){
    
    initLogger();
    
    // check arguments
    if (!checkParams(args)) {
      logger.trace("Exiting main method.");
      return;
    }

    String filename = args[0];
    logger.info("Graph file: " + filename);

    // initialize graphlab core
    final Core core;
    try {
      logger.trace("Initializing GraphLab core ...");
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(Scheduler.SWEEP);
      core = new Core(config);
    } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      logger.trace("Exiting main method.");
      return;
    }

    // construct graph
    logger.trace("Constructing graph from " + filename + " ...");
    final DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph;
    try {
      graph = constructGraph(filename);
    } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      logger.trace("Exiting main method.");
      core.destroy();
      return;
    }
    
    // execute graph updates
    core.setGraph(graph);
    core.scheduleAll(new PageRankUpdater(graph, 1));
    logger.trace("Running graphlab ...");
    logger.info("Took " + core.start() + " seconds.");
    
    // print results
    printResults (graph);
    logger.info("Update count: " + core.lastUpdateCount());
    
    // done
    core.destroy();
    
    return;

  }
  
  /** Initialize logger and set logging levels */
  private static void initLogger(){
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.INFO);
    logger.setLevel(Level.INFO);
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

  private static DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>
    constructGraph(String filename)
    throws IOException {

    DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    normalize(graph); 
    
    return graph;

  }
  
  private static void normalize(DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph){
    
    for(ScalarVertex vertex : graph.vertexSet()) {
      double sum = 0;
      Collection<DefaultWeightedEdge> outEdges = graph.outgoingEdgesOf(vertex);
      // Sum up weight on out edges
      for(DefaultWeightedEdge edge : outEdges){
        sum += graph.getEdgeWeight(edge);
      }
      for(DefaultWeightedEdge edge : outEdges){
        graph.setEdgeWeight(edge, (graph.getEdgeWeight(edge)/sum));
      }
    }
    
  }
  
  private static void printResults (DirectedGraph<ScalarVertex, DefaultWeightedEdge> g){
      
    logger.info("----------------- Results -----------------");
    logger.info("ID : Rank");
    
    Collection<ScalarVertex> vertices = g.vertexSet();
    List<ScalarVertex> verticesList = new ArrayList<ScalarVertex>(vertices.size());
    for (ScalarVertex vertex : vertices){
      verticesList.add(vertex);
    }
    
    Collections.sort(verticesList, new Comparator<ScalarVertex>(){
      public int compare(ScalarVertex left, ScalarVertex right) {
        return Double.compare(right.value(), left.value());
      }
    });
    
    for (int i=0; i<Math.min(verticesList.size(), 5); i++){
      ScalarVertex vertex = verticesList.get(i);
      logger.info(vertex.id() + " : " + vertex.value());
    }
    
  }

  private static class PageRankUpdater extends Updater<ScalarVertex> {

    /** Global reset probability */
    public static final double RESET_PROB = 0.15;
    
    /** Global accuracy tolerance */
    public static final double ACCURACY = 1e-5;
    
    /** Priority */
    private double mPriority;
    
    private DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> mGraph;

    public PageRankUpdater(
        DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph,
        double priority) {
      if (null == graph) throw new NullPointerException ("graph must not be null.");
      mGraph = graph;
      mPriority = priority;
    }

    @Override
    public void update(Context context, ScalarVertex vertex) {
      
      // Compute weighted sum of neighbors
      double sum = 0;
      
      /* Iterate over edge_id_list and get source is slow in graph2 */
      for(DefaultWeightedEdge edge : mGraph.incomingEdgesOf(vertex)){
        double weight = mGraph.getEdgeWeight(edge);
        sum += weight * mGraph.getEdgeSource(edge).value();
      }
      
      // Add random reset probability
      double oldValue = vertex.value();
      vertex.setValue(RESET_PROB + (1-RESET_PROB)*sum);
      for(DefaultWeightedEdge edge : mGraph.outgoingEdgesOf(vertex)) {    
        double weight = mGraph.getEdgeWeight(edge);
        double residual = weight * Math.abs(vertex.value() - oldValue);
        // If the neighbor changed sufficiently add to scheduler.
        if(residual > ACCURACY) 
          context.schedule(
              mGraph.getEdgeTarget(edge),
              new PageRankUpdater(mGraph, residual)
          );
      }

    }
    
    @Override
    public double priority(){
      return mPriority;
    }
    
    @Override
    public void add(Updater<ScalarVertex> other){
      if (!(other instanceof PageRankUpdater))
        throw new IllegalStateException("incompatible updaters added.");
      mPriority += ((PageRankUpdater) other).priority();
    }

    @Override
    protected Updater<ScalarVertex> clone() {
      return new PageRankUpdater(mGraph, mPriority);
    }

  }
  
}
