package org.graphlab.demo;

import java.io.IOException;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.graphlab.Aggregator;
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
 * Shortest Path algorithm.
 * 
 * <p>
 * Demonstrates GraphLab Java. To run this class, provide a path to a tsv
 * file in the program arguments. Some examples are available in the
 * <tt>java_jni/test-graphs</tt> directory.
 * </p>
 * 
 * <p>To run this demo, use</p>
<pre>
  ant ShortestPath -Dfile=&lt;path to tsv file&gt;
</pre>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ShortestPath {

  private static final Logger logger = Logger.getLogger(ShortestPath.class);

  public static void main(String[] args) {

    initLogger();

    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];

    // initialize graphlab core
    final Core core;
    try {
      CoreConfiguration config = new CoreConfiguration();
      config.setScheduler(Scheduler.FIFO);
      core = new Core();
    } catch (CoreException e) {
      logger.fatal("Unable to initialize core. Terminating.", e);
      return;
    }

    // construct graph
    final DirectedGraph<ScalarVertex, DefaultWeightedEdge> graph;
    try {
      graph = constructGraph(filename);
    } catch (IOException e) {
      core.destroy();
      return;
    }

    for (ScalarVertex v : graph.vertexSet())
      v.setValue(Integer.MAX_VALUE);

    // start from root
    ScalarVertex root = graph.vertexSet().iterator().next();
    root.setValue(0);

    core.setGraph(graph);
    core.schedule(root, new ShortestPathUpdater(graph));
    logger.trace("Runtime: " + core.start() + " s");
    logger.trace("Update count: " + core.lastUpdateCount());

    core.addAggregator("agg", new ShortestPathAggregator(), 0);
    core.aggregateNow("agg");
    
    // destroy core
    core.destroy();

  }

  private static void initLogger() {
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.ALL);
    logger.setLevel(Level.ALL);
  }

  /**
   * Checks that required input parameters are available and valid. Prints
   * instructions if not all parameters were valid.
   */
  private static boolean checkParams(String[] args) {

    if (1 == args.length) return true;
    
    System.out.println("Please provide filename.");
    System.out.println("Usage: java -Djava.library.path=... "
        + ShortestPath.class.getCanonicalName() + " path/to/tsv/file");
    return false;
      
  }

  private static DirectedGraph<ScalarVertex, DefaultWeightedEdge>
    constructGraph(String filename)
    throws IOException {

    DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge> graph =
       new DefaultDirectedWeightedGraph<ScalarVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    GraphLoader.loadGraphFromTsvFile(graph, filename);
    return graph;

  }
  
  /**
   * Finds shortest path to a particular vertex
   * @author Jiunn Haur Lim
   */
  private static class ShortestPathUpdater
    extends Updater<ScalarVertex, DefaultWeightedEdge, ShortestPathUpdater> {

    private DirectedGraph<ScalarVertex, DefaultWeightedEdge> mGraph;

    public ShortestPathUpdater(DirectedGraph<ScalarVertex, DefaultWeightedEdge> g) {
      this.mGraph = g;
    }

    @Override
    public void update(Context context, ScalarVertex vertex) {
      
      // find shortest known distance into this node
      for (DefaultWeightedEdge edge : mGraph.incomingEdgesOf(vertex)) {
        ScalarVertex neighbor = mGraph.getEdgeSource(edge);
        double weight = mGraph.getEdgeWeight(edge);
        vertex.setValue(
            (float) Math.min(vertex.value(), neighbor.value() + weight)
        );
      }

      // reschedule any affected neighbors
      for (DefaultWeightedEdge edge : mGraph.outgoingEdgesOf(vertex)) {
        ScalarVertex neighbor = mGraph.getEdgeTarget(edge);
        double weight = mGraph.getEdgeWeight(edge);
        if (neighbor.value() > (vertex.value() + weight)) {
          context.schedule(neighbor, this);
        }
      }

    }

    @Override
    protected ShortestPathUpdater clone() {
      return new ShortestPathUpdater(mGraph);
    }

  }

  /**
   * Aggregates shortest path information: total distance, longest distance, and furthest vertex
   * @author Jiunn Haur Lim
   */
  private static class ShortestPathAggregator
    extends Aggregator<ScalarVertex, ShortestPathAggregator> {
  
    /** maximum distance */
    private double mMaxDist;
    
    /** total distance */
    private double mSumDist;
    
    /** furthest vertex */
    private int mFurthestVertex;
    
    @Override
    protected void exec(Context context, ScalarVertex vertex) {
      double dist = vertex.value();
      mSumDist += dist;
      if (dist > mMaxDist){
        mMaxDist = dist;
        mFurthestVertex = vertex.id();
      }
    }
    
    @Override
    protected void add (ShortestPathAggregator other) {
      // add distances
      mSumDist += other.mSumDist;
      // take max
      if (other.mMaxDist > mMaxDist){
        mMaxDist = other.mMaxDist;
        mFurthestVertex = other.mFurthestVertex;
      }
    }
    
    @Override
    protected void finalize(Context context) {
       // output results
       logger.info ("Total Distance:\t\t" + mSumDist + "\n" + 
                    "Longest Distance:\t" + mMaxDist + "\n" +
                    "Furthest Vertex:\t"  + mFurthestVertex + "\n");
    }

    @Override
    protected ShortestPathAggregator clone() {
      return new ShortestPathAggregator();
    }
    
  }
  
}
