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
    
    ///////// DUMB VISUAL VM HACK ///////////
    // Scanner sc = new Scanner(System.in);
    // sc.next();
    ///////// DUMB VISUAL VM HACK ///////////
    
    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];

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
    final DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph;
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
    core.scheduleAll(new PageRankUpdater(graph, PageRankUpdater.RESET_PROB));
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
    logger.setLevel(Level.ALL);
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
          + PageRank.class.getCanonicalName() + " path/to/tsv/file");
      return false;
    }

    return true;

  }

  private static DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge>
    constructGraph(String filename)
    throws IOException {

    DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    // GraphLoader.loadGraphFromAdjFile(graph, PageRankVertex.class, filename);
    GraphLoader.loadGraphFromTsvFile(graph, PageRankVertex.class, filename);
    normalize(graph); 
    
    return graph;

  }
  
  private static void normalize(DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph){
    
    for(PageRankVertex vertex : graph.vertexSet()) {
      double sum = 0;
      Collection<DefaultWeightedEdge> outEdges = graph.outgoingEdgesOf(vertex);
      // sum up weight on out edges
      for(DefaultWeightedEdge edge : outEdges){
        sum += graph.getEdgeWeight(edge);
      }
      for(DefaultWeightedEdge edge : outEdges){
        graph.setEdgeWeight(edge, (graph.getEdgeWeight(edge)/sum));
      }
    }
    
  }
  
  private static void printResults (DirectedGraph<PageRankVertex, DefaultWeightedEdge> g){
      
    logger.info("----------------- Results -----------------");
    logger.info("ID : Rank");
    
    Collection<PageRankVertex> vertices = g.vertexSet();
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

  private static class PageRankUpdater extends Updater<PageRankVertex, DefaultWeightedEdge, PageRankUpdater> {

    /** Global reset probability */
    public static final double RESET_PROB = 0.15;
    
    /** Global accuracy tolerance */
    public static final double ACCURACY = 1e-5;
    
    /** Priority */
    private double mPriority;
    
    private DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> mGraph;

    public PageRankUpdater(
        DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph,
        double priority) {
      if (null == graph) throw new NullPointerException ("graph must not be null.");
      mGraph = graph;
      mPriority = priority;
    }

    @Override
    public void update(Context context, PageRankVertex vertex) {
      
      vertex.mNUpdates++;
      
      // compute weighted sum of neighbors
      double sum = 0;
      
      /* Iterate over edge_id_list and get source is slow in graph2 */
      for(DefaultWeightedEdge edge : mGraph.incomingEdgesOf(vertex)){
        double weight = mGraph.getEdgeWeight(edge);
        sum += weight * mGraph.getEdgeSource(edge).value();
      }
      
      // add random reset probability
      vertex.mOldValue = vertex.value();
      vertex.setValue(RESET_PROB + (1-RESET_PROB)*sum);
      for(DefaultWeightedEdge edge : mGraph.outgoingEdgesOf(vertex)) {    
        double weight = mGraph.getEdgeWeight(edge);
        double residual = weight * Math.abs(vertex.value() - vertex.mOldValue);
        // if the neighbor changed sufficiently add to scheduler.
        if(residual > ACCURACY) 
          context.schedule(mGraph.getEdgeTarget(edge),
              new PageRankUpdater(mGraph, residual));
      }

    }
    
    @Override
    public double priority(){
      return mPriority;
    }
    
    @Override
    public void add(PageRankUpdater other){
      mPriority += other.priority();
    }

    @Override
    protected PageRankUpdater clone() {
      return new PageRankUpdater(mGraph, mPriority);
    }

  }
  
}
