package org.graphlab.toolkits.matrix.als;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import no.uib.cipr.matrix.Matrices;
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
import org.graphlab.Core;
import org.graphlab.Core.CoreException;
import org.graphlab.CoreConfiguration;
import org.graphlab.Scheduler;
import org.graphlab.toolkits.matrix.MatrixLoader;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Matrix factorization w. alternating least squares.
 * Adapted from Gonzalez's original ALS code.
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a
 *      href="http://code.google.com/p/graphlabapi/source/browse/toolkits/matrix_factorization/als.cpp?name=v2">GraphLab
 *      ALS</a>
 */
public class Als {

  /** Maximum possible rating */
  private static double UPPER_BOUND = Double.MAX_VALUE;

  /** Minimum possible rating */
  private static double LOWER_BOUND = Double.MIN_VALUE;

  protected static final Logger logger = Logger.getLogger(Als.class);

  /**
   * Executes ALS matrix factorization on input matrix.
   * 
   * @param args
   *          Command line arguments. Use `java ... -h` to see what they are.
   * @throws IOException
   * @throws CoreException
   * @throws IllegalAccessException
   * @throws InstantiationException
   */
  public static void main(String[] args)
    throws IOException,
           CoreException,
           InstantiationException,
           IllegalAccessException {

    // parse command line ------------------------------------------------------
    CommandLine cmd = parseCommandLine(args);
    String filename = cmd.getOptionValue("train-data");
    if (cmd.hasOption("upper-bound"))
      UPPER_BOUND = Double.parseDouble(cmd.getOptionValue("upper-bound"));
    if (cmd.hasOption("lower-bound"))
      LOWER_BOUND = Double.parseDouble(cmd.getOptionValue("lower-bound"));

    initLogger();

    // construct graph ---------------------------------------------------------
    Map<Integer, AlsVertex> vertices = new HashMap<Integer, AlsVertex>();
    AlsGraph graph = new AlsGraph(DefaultWeightedEdge.class);
    MatrixLoader.loadGraphFromMM(graph, vertices, AlsVertex.class, filename);
    
    // initialize random latent factors ----------------------------------------
    randomLatentFactors(graph, AlsUpdater.NLATENT);

    // init graphlab core ------------------------------------------------------
    final CoreConfiguration config = new CoreConfiguration();
    config.setScheduler(Scheduler.SWEEP);
    final Core core = new Core(config);

    // schedule and run --------------------------------------------------------
    core.setGraph(graph);
    core.scheduleAll(new AlsUpdater(graph));
    logger.info("GraphLab engine stopped. Took " + core.start() + " seconds.");

    // aggregate statistics ----------------------------------------------------
    core.addAggregator("agg", new AlsAggregator(graph), 0);
    core.aggregateNow("agg");

    // done --------------------------------------------------------------------
    core.destroy();

    // validate ----------------------------------------------------------------
    if (cmd.hasOption("test-data"))
      validate(vertices, cmd.getOptionValue("test-data"));

  }

  /**
   * Parses command line arguments. Prints help message and exits on
   * ParseException.
   * 
   * @param args      command line arguments
   * @return parser
   */
  @SuppressWarnings("static-access")
  private static CommandLine parseCommandLine(String[] args) {

    CommandLineParser parser = new PosixParser();
    Options options = new Options();
    CommandLine cmd = null;

    try {

      Option trainData = OptionBuilder.withArgName("train-data")
          .withDescription("File containing training matrix")
          .withLongOpt("train-data").withType("").isRequired().hasArg()
          .create();
      Option testData = OptionBuilder.withArgName("test-data")
          .withDescription("File containing test matrix")
          .withLongOpt("test-data").withType("").hasArg().create();
      Option upperBound = OptionBuilder.withArgName("upper-bound")
          .withDescription("Maximum rating").withLongOpt("upper-bound")
          .withType(0.0).hasArg().create();
      Option lowerBound = OptionBuilder.withArgName("lower-bound")
          .withDescription("Minimum rating").withLongOpt("lower-bound")
          .withType(0.0).hasArg().create();

      options.addOption(trainData);
      options.addOption(testData);
      options.addOption(upperBound);
      options.addOption(lowerBound);

      cmd = parser.parse(options, args);

    } catch (ParseException e) {
      // automatically generate the help statement
      new HelpFormatter().printHelp(Als.class.getCanonicalName(), options);
      System.exit(2);
    }

    return cmd;

  }

  /**
   * Initializes vertices with random latent factors between 0 and 1.
   * 
   * @param graph
   *          graph to populate
   * @param nlatent
   *          number of latent factors
   */
  private static void randomLatentFactors(AlsGraph graph, int nlatent) {
    for (AlsVertex vertex : graph.vertexSet()) {
      vertex.setVector(Matrices.random(nlatent));
    }
  }
  
  /** Initialize logger and set logging levels */
  private static void initLogger() {
    BasicConfigurator.configure();
    Logger.getLogger(Core.class).setLevel(Level.INFO);
    logger.setLevel(Level.ALL);
  }

  /**
   * Validates learned model against training data. Prints test error.
   * 
   * <p>
   * After matrix factorization is completed, every vertex contains a vector
   * of latent factors. By taking the dot product of a user vector and an 
   * item vector, we can predict the rating that the user might give to the
   * item.
   * </p>
   * 
   * @param vertices
   *          Map containing ALS vertices; each vertex contains a vector of
   *          latent factors.
   * @param filename
   *          Path to test file.
   * @throws IOException
   *          if it is unable to parse the test data
   */
  private static void validate(Map<Integer, AlsVertex> vertices, String filename)
      throws IOException {

    // read matrix meta-data
    MatrixVectorReader reader = new MatrixVectorReader(new FileReader(filename));
    MatrixSize size = reader.readMatrixSize(reader.readMatrixInfo());

    int[] row = new int[1];
    int[] col = new int[1];
    double[] data = new double[1];
    
    // calculate mean
    double testMean = 0;
    for (int i = 0; i < size.numEntries(); i++){
      reader.readCoordinate(row, col, data);
      testMean += data[0];
    }
    testMean /= size.numEntries();

    double sumSquaredErrors = 0;
    reader = new MatrixVectorReader(new FileReader(filename));
    size = reader.readMatrixSize(reader.readMatrixInfo());

    // iterate through file entries and test
    for (int i = 0; i < size.numEntries(); i++) {

      reader.readCoordinate(row, col, data);

      int sourceID = row[0];
      int targetID = size.numRows() + col[0];
      double expected = data[0];

      AlsVertex source = vertices.get(sourceID);
      AlsVertex target = vertices.get(targetID);
      if (null == source)
        throw new IOException("Vertex " + sourceID + " not found.");
      if (null == target)
        throw new IOException("Vertex " + targetID + " not found.");

      double predicted = truncate(source.vector().dot(target.vector()));
      double error = predicted - expected;
      sumSquaredErrors += error * error;

    }

    double rmse = Math.sqrt(sumSquaredErrors / size.numEntries());
    logger.info("-------- Test Results --------");
    logger.info("RMSE: " + rmse);

  }
  
  /**
   * Constraints <tt>value</tt> within {@link #LOWER_BOUND} and
   * {@link #UPPER_BOUND}.
   * @param value
   * @return
   *  if LOWER < value < UPPER, returns value.<br/>
   *  if value > UPPER, returns UPPER.<br/>
   *  if value < LOWER, returns LOWER.<br/.
   */
  protected static double truncate(double value){
    value = Math.min(value, UPPER_BOUND);
    value = Math.max(value, LOWER_BOUND);
    return value;
  }

}
