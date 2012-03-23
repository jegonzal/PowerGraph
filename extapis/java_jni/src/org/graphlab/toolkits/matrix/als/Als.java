package org.graphlab.toolkits.matrix.als;

import java.io.IOException;
import java.util.Set;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
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

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleCholeskyDecomposition;


/**
 * Matrix factorization w. alternating least squares.
 * TODO allow edge to store observation and weights
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * Adapted from Danny Bickson's original ALS code
 */
public class Als {
  
  /** Number of latent factors */
  // TODO
  // private static final int NLATENT = 20;
  private static final int NLATENT = 2;
  
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
    MatrixLoader.loadGraph(graph, AlsVertex.class, filename);
    randomLatentFactors(graph, NLATENT);
    
    // init graphlab core
    final CoreConfiguration config = new CoreConfiguration();
    config.setScheduler(Scheduler.SWEEP);
    config.setNCpus(2); // TODO
    final Core core = new Core(config);
    
    // schedule and run
    core.setGraph(graph);
    core.scheduleAll(new AlsUpdater(graph, 10000));
    logger.info("GraphLab engine stopped. Took " + core.start() + " seconds.");
    
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
    
    DoubleFactory1D factory = DoubleFactory1D.dense;
    for (VectorVertex vertex : graph.vertexSet())
      vertex.setVector(factory.random(nlatent));
    
  }
  
  /** Initialize logger and set logging levels */
  private static void initLogger(){
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.INFO);
    logger.setLevel(Level.ALL);
  }
  
  /**
   * Alternating Least Squares updater. TODO
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
      
      logger.debug("Updating vertex " + vertex.id());
      
      vertex.mSquaredError = vertex.mResidual = 0;
      ++vertex.mNUpdates;
      
      // if there are no neighbors just return
      Set<DefaultWeightedEdge> edges = mGraph.edgesOf(vertex);
      if(edges.isEmpty()) return;
      
      // get the number of latent dimensions
      DoubleMatrix2D XtX = DoubleFactory2D.dense.make(NLATENT, NLATENT);
      DoubleMatrix1D Xty = DoubleFactory1D.dense.make(NLATENT);
        
      // Compute X'X and X'y (weighted) -----------------------------------------
      for (final DefaultWeightedEdge edge : edges) {
        // get neighbor
        DoubleMatrix1D neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex).vector();
        double weight = mGraph.getEdgeWeight(edge);
        // update the X'X and X'y
        for(int i = 0; i < NLATENT; ++i) {
          double increment = neighbor.getQuick(i)*weight; // weight is the y value
          Xty.setQuick(i, Xty.getQuick(i) + increment);
          for(int j = i; j < NLATENT; ++j){
            increment = neighbor.getQuick(i)*neighbor.getQuick(j);
            XtX.setQuick(j, i, XtX.getQuick(j, i) + increment);
          }
        }
      }
      
      for(int i = 0; i < NLATENT; ++i)
        for(int j = i+1; j < NLATENT; ++j)
          XtX.setQuick(i, j, XtX.getQuick(j, i));

      for(int i = 0; i < NLATENT; ++i)
        XtX.setQuick(i, i, XtX.getQuick(i,i) + (LAMBDA) * edges.size());
      
      logger.debug(XtX);
      logger.debug(Xty);
      
      // solve the least squares problem using PColt ----------------------------
      final DoubleMatrix1D oldLatent = vertex.vector();
      new DenseDoubleCholeskyDecomposition(XtX).solve(Xty); // in-place
      vertex.setVector(Xty);

      logger.debug(Xty);
      
      // compute the residual change in the latent factors-----------------------
      vertex.mResidual = 0;
      for(int i = 0; i < NLATENT; ++i)
        vertex.mResidual += Math.abs(oldLatent.getQuick(i) - vertex.vector().getQuick(i));
      vertex.mResidual /= NLATENT;
      
      // update the rmse and reschedule neighbors -------------------------------
      logger.debug("Residual: " + vertex.mResidual);
      for (final DefaultWeightedEdge edge : edges) {
        // get the neighbor id
        AlsVertex neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex);
        final double pred = vertex.vector().zDotProduct(neighbor.vector());
        final double error = Math.abs(mGraph.getEdgeWeight(edge) - pred);
        vertex.mSquaredError += error*error;
        // reschedule neighbors ------------------------------------------------
        logger.debug("\t" + neighbor.id() + ": " + error);
        if( error > TOLERANCE && vertex.mResidual > TOLERANCE) 
          context.schedule(neighbor, new AlsUpdater(mGraph, error * vertex.mResidual));
      }
      
      logger.debug("scheduled.");
      
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

}
