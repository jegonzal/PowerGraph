package org.graphlab.demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Scanner;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Context;
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.Updater;
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
public class PageRankFactorized {

  private static final Logger logger = Logger.getLogger(PageRankFactorized.class);
  
  public static void main (String[] args){
    
    initLogger();
    
    ///////// DUMB VISUAL VM HACK ///////////
    Scanner sc = new Scanner(System.in);
    sc.next();
    ///////// DUMB VISUAL VM HACK ///////////
    
    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];

    // initialize graphlab core
    final Core core;
    try {
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(Scheduler.SWEEP);
      core = new Core(config);
    } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      return;
    }

    // construct graph
    final DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph;
    try {
      graph = constructGraph(filename);
    } catch (IOException e) {
      logger.fatal("Unable to construct graph. Terminating.", e);
      core.destroy();
      return;
    }
    
    // execute graph updates
    core.setGraph(graph);
    core.scheduleAll(new PageRankUpdater(graph, PageRankUpdater.RESET_PROB));
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

  /**
   * Constructs graph and normalizes weights
   * @param filename
   * @return graph
   * @throws IOException
   */
  private static DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge>
    constructGraph(String filename)
    throws IOException {

    DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph
      = new DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
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
      // set normalized weight
      for(DefaultWeightedEdge edge : outEdges){
        graph.setEdgeWeight(edge, (graph.getEdgeWeight(edge)/sum));
      }
    }
    
  }
  
  /**
   * Prints results for top 5 pages.
   * @param graph
   */
  private static void printResults (DirectedGraph<PageRankVertex, DefaultWeightedEdge> graph){
      
    logger.info("----------------- Results -----------------");
    logger.info("ID : Rank");
    
    Collection<PageRankVertex> vertices = graph.vertexSet();
    List<PageRankVertex> verticesList = new ArrayList<PageRankVertex>(vertices.size());
    for (PageRankVertex vertex : vertices){
      verticesList.add(vertex);
    }
    
    Collections.sort(verticesList, new Comparator<PageRankVertex>(){
      public int compare(PageRankVertex left, PageRankVertex right) {
        return Double.compare(right.value(), left.value());
      }
    });
    
    // print top 5 pages
    for (int i=0; i<Math.min(verticesList.size(), 5); i++){
      PageRankVertex vertex = verticesList.get(i);
      logger.info(vertex.id() + " : " + vertex.value());
    }
    
  }

  private static class PageRankUpdater
    extends Updater<PageRankVertex, DefaultWeightedEdge, PageRankUpdater> {

    /** Global reset probability */
    public static final double RESET_PROB = 0.15;
    
    /** Global accuracy tolerance */
    public static final double ACCURACY = 1e-5;
    
    /** Accumulated PageRank */
    private double mAccum;
    
    private DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> mGraph;

    public PageRankUpdater(
        DefaultDirectedWeightedGraph<PageRankVertex, DefaultWeightedEdge> graph,
        double accum) {
      mGraph = graph;
      mAccum = accum;
    }
    
    @Override
    public double priority(){
      return Math.abs(mAccum);
    }
    
    @Override
    public void add(PageRankUpdater other){
      mAccum += other.priority();
    }

    @Override
    protected PageRankUpdater clone() {
      return new PageRankUpdater(mGraph, mAccum);
      // return this; // this is really not the problem
    }
    
    @Override
    protected boolean isFactorizable(){
      return true;
    }
    
    @Override
    protected int gatherConsistency(){
      return Updater.EDGE_CONSISTENCY;
    }
    
    @Override
    protected int scatterConsistency(){
      return Updater.NULL_CONSISTENCY;
    }
    
    @Override
    protected int gatherEdges(){
      return Updater.IN_EDGES;
    }
    
    @Override
    protected int scatterEdges() {
      return (Math.abs(mAccum) > ACCURACY) ?
          Updater.OUT_EDGES : Updater.NO_EDGES;
    }
    
    // reset the accumulator before running the gather
    @Override
    protected void initGather(){
      mAccum = 0;
    }
    
    // run the gather operation over all in edges
    @Override
    protected void gather(DefaultWeightedEdge edge) {
      mAccum +=
        mGraph.getEdgeSource(edge).value() *
        mGraph.getEdgeWeight(edge);
    }
    
    @Override
    protected void merge(PageRankUpdater updater){
      // accumulate page ranks
      mAccum += ((PageRankUpdater) updater).mAccum;
    }
    
    @Override
    protected void apply(PageRankVertex vertex) {
      // update the center vertex
      vertex.mNUpdates++;
      vertex.setValue(RESET_PROB + (1 - RESET_PROB) * mAccum);
      mAccum = vertex.value() - vertex.mOldValue;
      if(Math.abs(mAccum) > ACCURACY || vertex.mNUpdates == 1) {
        vertex.mOldValue = vertex.value();    
      }
    }

    @Override
    protected void scatter(Context context, DefaultWeightedEdge edge) {
      // reschedule neighbors 
      double delta = mAccum * mGraph.getEdgeWeight(edge) * (1 - RESET_PROB);
      context.schedule(mGraph.getEdgeTarget(edge),
          new PageRankUpdater(mGraph, delta));
    }

  }
  
}
