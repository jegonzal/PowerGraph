package org.graphlab.toolkits.matrix.als;

import java.io.IOException;
import java.util.Set;

import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

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
import org.graphlab.toolkits.matrix.MatrixLoader;
import org.graphlab.toolkits.matrix.VectorVertex;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;


/**
 * Matrix factorization w. alternating least squares.
 * TODO allow edge to store observation and weights
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * Adapted from Danny Bickson's original ALS code
 */
public class Als {
  
  /** Number of latent factors */
  private static final int NLATENT = 20;
  // private static final int NLATENT = 2;
  
  private static final double LAMBDA = 0.065;
  private static final double TOLERANCE = 1e-2;
  
  private static final Logger logger = Logger.getLogger(Als.class);
  
  /**
   * Executes ALS matrix factorization on input file.
   * @param args
   *          command line arguments
   * @throws IOException 
   * @throws CoreException 
   * @throws IllegalAccessException 
   * @throws InstantiationException 
   */
  public static void main(String[] args)
    throws IOException, CoreException, InstantiationException, IllegalAccessException {
    
    // check arguments
    if (!checkParams(args)) return;
    String filename = args[0];
    
    initLogger();
    
    // construct graph
    SimpleWeightedGraph<AlsVertex, DefaultWeightedEdge>
      graph = new SimpleWeightedGraph<AlsVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    MatrixLoader.loadGraphFromMM(graph, AlsVertex.class, filename);
    randomLatentFactors(graph, NLATENT);
    
    // init graphlab core
    final CoreConfiguration config = new CoreConfiguration();
    config.setScheduler(Scheduler.SWEEP);
    config.setNCpus(8);
    final Core core = new Core(config);
    
    // schedule and run
    core.setGraph(graph);
    core.scheduleAll(new AlsUpdater(graph, 10000));
    logger.info("GraphLab engine stopped. Took " + core.start() + " seconds.");
    
    // aggregate statistics
    core.addAggregator("agg", new AlsAggregator(graph), 0);
    core.aggregateNow("agg");
    
    core.destroy();
    
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
          + Als.class.getCanonicalName() + " path/to/matrix/file");
      return false;
    }

    return true;

  }
  
  /**
   * Initialize vertices with random latent factors.
   * @param graph
   * @param nlatent
   *          number of latent factors
   */
  private static void
    randomLatentFactors(Graph<AlsVertex, ?> graph, int nlatent){
    
    for (VectorVertex vertex : graph.vertexSet()){
      double[] values = new double[nlatent];
      for (int i=0; i<values.length; i++) values[i] = Math.random();
      vertex.setVector(new DenseVector(values));
    }
    
  }
  
  /** Initialize logger and set logging levels */
  private static void initLogger(){
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.INFO);
    logger.setLevel(Level.ALL);
  }
  
  /**
   * Alternating Least Squares updater.
   * 
   * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
   */
  private static class AlsUpdater
    extends Updater<AlsVertex, DefaultWeightedEdge, AlsUpdater> {

    private double mError;
    private WeightedGraph<AlsVertex, DefaultWeightedEdge> mGraph;
    
    public AlsUpdater(WeightedGraph<AlsVertex, DefaultWeightedEdge> graph, double error){
      mError = error;
      mGraph = graph;
    }
    
    @Override
    protected AlsUpdater clone() {
      return new AlsUpdater(mGraph, mError);
    }
    
    @Override
    public void update(Context context, AlsVertex vertex) {
      
      vertex.mSquaredError = vertex.mResidual = 0;
      ++vertex.mNUpdates;
      
      // if there are no neighbors just return
      Set<DefaultWeightedEdge> edges = mGraph.edgesOf(vertex);
      if(edges.isEmpty()) return;
      
      // get the number of latent dimensions
      DenseMatrix XtX = new DenseMatrix(NLATENT, NLATENT);
      DenseVector Xty = new DenseVector(NLATENT);
        
      // Compute X'X and X'y (weighted) -----------------------------------------
      for (final DefaultWeightedEdge edge : edges) {
        // get neighbor
        DenseVector neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex).vector();
        double weight = mGraph.getEdgeWeight(edge);
        // update the X'X and X'y
        for(int i = 0; i < NLATENT; ++i) {
          double increment = neighbor.get(i)*weight; // weight is the y value
          Xty.set(i, Xty.get(i) + increment);
          for(int j = i; j < NLATENT; ++j){
            increment = neighbor.get(i)*neighbor.get(j);
            XtX.set(j, i, XtX.get(j, i) + increment);
          }
        }
      }
      
      for(int i = 0; i < NLATENT; ++i)
        for(int j = i+1; j < NLATENT; ++j)
          XtX.set(i, j, XtX.get(j, i));

      for(int i = 0; i < NLATENT; ++i)
        XtX.set(i, i, XtX.get(i,i) + (LAMBDA) * edges.size());
      
      // solve the least squares problem using PColt ----------------------------
      final DenseVector oldLatent = vertex.vector();
      DenseMatrix newLatent = DenseCholesky.factorize(XtX).solve(new DenseMatrix(Xty));
      vertex.setVector(new DenseVector(newLatent.getData()));
      
      // compute the residual change in the latent factors-----------------------
      vertex.mResidual = 0;
      for(int i = 0; i < NLATENT; ++i)
        vertex.mResidual += Math.abs(oldLatent.get(i) - vertex.vector().get(i));
      vertex.mResidual /= NLATENT;
      
      // update the rmse and reschedule neighbors -------------------------------
      for (final DefaultWeightedEdge edge : edges) {
        // get the neighbor id
        AlsVertex neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex);
        final double pred = vertex.vector().dot(neighbor.vector());
        final double error = Math.abs(mGraph.getEdgeWeight(edge) - pred);
        vertex.mSquaredError += error*error;
        // reschedule neighbors ------------------------------------------------
        if( error > TOLERANCE && vertex.mResidual > TOLERANCE) 
          context.schedule(neighbor, new AlsUpdater(mGraph, error * vertex.mResidual));
      }
      
    } // end of operator()
    
    @Override
    public void add(AlsUpdater other){
      mError += other.mError;
    }
    
    @Override
    public double priority(){
      return mError;
    }
    
  }
  
  /**
   * Aggregate results of ALS
   * @author Jiunn Haur Lim
   */
  private static class AlsAggregator
    extends Aggregator<AlsVertex, AlsAggregator>{

    private double mSumSquaredErrors = 0;
    private double mMaxRootMeanSquaredError = 0;
    private WeightedGraph<AlsVertex, DefaultWeightedEdge> mGraph;
    
    private AlsAggregator(WeightedGraph<AlsVertex, DefaultWeightedEdge> graph){
      mGraph = graph;
    }
    
    @Override
    protected void exec(Context context, AlsVertex vertex) {
      
      int numEdges = mGraph.edgesOf(vertex).size();
      if (0 == numEdges) return;
      
      mSumSquaredErrors += vertex.mSquaredError;
      mMaxRootMeanSquaredError =
          Math.max(mMaxRootMeanSquaredError, Math.sqrt(vertex.mSquaredError/numEdges));
      
    }

    @Override
    protected void add(AlsAggregator other) {
      mSumSquaredErrors += other.mSumSquaredErrors;
      mMaxRootMeanSquaredError =
        Math.max(mMaxRootMeanSquaredError, other.mMaxRootMeanSquaredError);
    }

    @Override
    protected void finalize(Context context) {
      logger.info("Average RMS: " +
          Math.sqrt(mSumSquaredErrors/(2*mGraph.edgeSet().size())));
      logger.info("Max RMS: " + mMaxRootMeanSquaredError);
    }

    @Override
    protected AlsAggregator clone() {
      AlsAggregator agg = new AlsAggregator(mGraph);
      agg.mMaxRootMeanSquaredError = this.mMaxRootMeanSquaredError;
      agg.mSumSquaredErrors = this.mSumSquaredErrors;
      return agg;
    }
    
  }

}
