package org.graphlab.toolkits.matrix.als;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
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
 * Adapted from Gonzalez's original ALS code.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see
 * <a href="http://code.google.com/p/graphlabapi/source/browse/toolkits/matrix_factorization/als.cpp?name=v2">GraphLab ALS</a>
 */
public class Als {
  
  /** Number of latent factors */
  private static final int NLATENT = 20;
  
  /** Regularization parameter */
  private static final double LAMBDA = 0.065;
  
  /** Convergence tolerance */
  private static final double TOLERANCE = 1e-2;
  
  /** Maximum possible rating */
  private static double UPPER_BOUND = Double.MAX_VALUE;
  
  /** Minimum possible rating */
  private static double LOWER_BOUND = Double.MIN_VALUE;
  
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
    
    // parse command line
    CommandLine cmd = parseCommandLine(args);
    String filename = cmd.getOptionValue("train-data");
    
    if (cmd.hasOption("upper-bound")){
      UPPER_BOUND = Double.parseDouble(cmd.getOptionValue("upper-bound"));
    }
    
    if (cmd.hasOption("lower-bound")){
      LOWER_BOUND = Double.parseDouble(cmd.getOptionValue("lower-bound"));
    }
    
    initLogger();
    
    // construct graph
    Map<Integer, AlsVertex> vertices = new HashMap<Integer, AlsVertex>();
    SimpleWeightedGraph<AlsVertex, DefaultWeightedEdge>
      graph = new SimpleWeightedGraph<AlsVertex, DefaultWeightedEdge>(DefaultWeightedEdge.class);
    MatrixLoader.loadGraphFromMM(graph, vertices, AlsVertex.class, filename);
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
    
    // done
    core.destroy();
    
    // validate
    if (cmd.hasOption("test-data")){
      validate(vertices, cmd.getOptionValue("test-data"));
    }
    
  }
  
  /**
   * Parses command line arguments. Prints help message and exits on
   * ParseException.
   * @param args
   * @return parser
   */
  @SuppressWarnings("static-access")
  private static CommandLine parseCommandLine(String[] args){
    
    CommandLineParser parser = new PosixParser();
    Options options = new Options();
    CommandLine cmd = null;
    
    try {
      
      Option trainData = OptionBuilder
                           .withArgName("train-data")
                           .withDescription("File containing training matrix")
                           .withLongOpt("train-data")
                           .withType("")
                           .isRequired()
                           .hasArg()
                           .create();
      Option testData = OptionBuilder
                          .withArgName("test-data")
                          .withDescription("File containing test matrix")
                          .withLongOpt("test-data")
                          .withType("")
                          .hasArg()
                          .create();
      Option upperBound = OptionBuilder
                            .withArgName("upper-bound")
                            .withDescription("Maximum rating")
                            .withLongOpt("upper-bound")
                            .withType(0.0)
                            .hasArg()
                            .create();
      Option lowerBound = OptionBuilder
                            .withArgName("lower-bound")
                            .withDescription("Minimum rating")
                            .withLongOpt("lower-bound")
                            .withType(0.0)
                            .hasArg()
                            .create();
      
      options.addOption(trainData);
      options.addOption(testData);
      options.addOption(upperBound);
      options.addOption(lowerBound);
      
      cmd = parser.parse(options, args);
    
    }catch (ParseException e){
      // automatically generate the help statement
      new HelpFormatter().printHelp( Als.class.getCanonicalName(), options );
      System.exit(2);
    }
    
    return cmd;
    
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
   * Validates learned model against training data. Prints test error.
   * @param vertices
   *          Map containing ALS vertices; each vertex contains a vector of latent factors.
   * @param filename
   *          Path to test file.
   * @throws IOException
   */
  private static void
    validate(Map<Integer, AlsVertex> vertices, String filename)
    throws IOException{
    
    // read matrix meta-data
    MatrixVectorReader reader = new MatrixVectorReader(new FileReader(filename));
    MatrixSize size = reader.readMatrixSize(reader.readMatrixInfo());
    
    int[] row = new int[1];
    int[] col = new int[1];
    double[] data = new double[1];
    double sumSquaredErrors = 0;
    
    // iterate through file entries and test
    for (int i=0; i<size.numEntries(); i++){
      
      reader.readCoordinate(row, col, data);
      
      int sourceID = row[0];
      int targetID = size.numRows() + col[0];
      double expected = data[0];
      
      AlsVertex source = vertices.get(sourceID);
      AlsVertex target = vertices.get(targetID);
      if (null == source)
        throw new IOException ("Vertex " + sourceID + " not found.");
      if (null == target)
        throw new IOException ("Vertex " + targetID + " not found.");
      
      double error = source.vector().dot(target.vector()) - expected;
      sumSquaredErrors = error * error;
      
    }
    
    double rms = sumSquaredErrors/size.numEntries();
    logger.info("-------- Test Results --------");
    logger.info("RMS: " + rms);
    
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
        double pred = vertex.vector().dot(neighbor.vector());
        
        // truncate
        pred = Math.min(pred, UPPER_BOUND);
        pred = Math.max(pred, LOWER_BOUND);
        
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
      logger.info("-------- Train Results --------");
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
