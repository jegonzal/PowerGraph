package org.graphlab.demo;

import java.io.IOException;
import java.util.Collection;

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
    final SparseGraph<ScalarVertex, ScalarEdge> graph;
    try {
      graph = constructGraph(filename);
    } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      logger.trace("Exiting main method.");
      c.destroy();
      return;
    }
    
    // execute graph updates
    c.setGraph(graph);
    c.scheduleAll(new PageRankUpdater(c, graph, 1));
    logger.trace("Running graphlab ...");
    logger.info("Took " + c.start() + " seconds.");
    
    // print results
    printResults (graph);
    
    // done
    c.destroy();
    
    return;

  }
  
  /** Initialize logger and set logging levels */
  private static void initLogger(){
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
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
          + ShortestPath.class.getCanonicalName() + " path/to/tsv/file");
      return false;
    }

    return true;

  }

  private static SparseGraph<ScalarVertex, ScalarEdge> constructGraph(
      String filename) throws IOException {

    SparseGraph<ScalarVertex, ScalarEdge> graph = new SparseGraph<ScalarVertex, ScalarEdge>();
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    normalize(graph); 
    
    return graph;

  }
  
  private static void normalize(SparseGraph<ScalarVertex, ScalarEdge> graph){
    for(ScalarVertex vertex : graph.vertices()) {
      double sum = 0;
      Collection<ScalarEdge> outEdges = graph.outgoingEdges(vertex.id());
      // Sum up weight on out edges
      for(ScalarEdge edge : outEdges){
        sum += edge.weight();
      }
      for(ScalarEdge edge : outEdges){
        edge.setWeight(edge.weight()/sum);
      }
    }
  }
  
  private static void printResults (SparseGraph<ScalarVertex, ScalarEdge> g){
      
    logger.info("----------------- Results -----------------");
    logger.info("ID : Color");
    for (ScalarVertex v : g.vertices()){
      logger.info(v.id() + " : " + (int) v.value());
    }
    
  }

  private static class PageRankUpdater extends Updater {

    /** Global reset probability */
    public static final double RESET_PROB = 0.15;
    
    /** Global accuracy tolerance */
    public static final double ACCURACY = 1e-5;
    
    /** Priority */
    private double mPriority;
    
    private SparseGraph<ScalarVertex, ScalarEdge> mGraph;
    private Core<?> mCore;

    public PageRankUpdater(Core<?> core, SparseGraph<ScalarVertex, ScalarEdge> graph, double priority) {
      super(core);
      if (null == graph) throw new NullPointerException ("graph must not be null.");
      mCore = core;
      mGraph = graph;
      mPriority = priority;
    }

    @Override
    public void update(Context context, int vertexId) {

      ScalarVertex vertex = mGraph.getVertex(vertexId);
      
      // Compute weighted sum of neighbors
      double sum = 0;
      
      /* Iterate over edge_id_list and get source is slow in graph2 */
      for(ScalarEdge edge : mGraph.incomingEdges(vertexId)){ 
        sum += edge.weight() * 
          mGraph.getVertex(edge.source()).value();
      }
      
      // Add random reset probability
      double oldValue = vertex.value();
      vertex.setValue(RESET_PROB + (1-RESET_PROB)*sum);
      for(ScalarEdge edge : mGraph.outgoingEdges(vertexId)) {    
        double residual = edge.weight() * Math.abs(vertex.value()- oldValue);
        // If the neighbor changed sufficiently add to scheduler.
        if(residual > ACCURACY) 
          context.schedule(edge.target(), new PageRankUpdater(mCore, mGraph, residual));
      }

    }
    
    @Override
    public double priority(){
      return mPriority;
    }
    
    @Override
    public void add(Updater other){
      if (!(other instanceof PageRankUpdater))
        throw new IllegalStateException("incompatible updaters added.");
      mPriority += other.priority();
    }

  }
  
}
