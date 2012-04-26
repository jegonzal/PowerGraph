package org.graphlab.toolkits.matrix.als;

import java.util.Set;

import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.graphlab.Context;
import org.graphlab.Updater;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Alternating Least Squares updater.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class AlsUpdater
extends Updater<AlsVertex, DefaultWeightedEdge, AlsUpdater> {
  
  /** Regularization parameter */
  private static final double LAMBDA = 0.065;

  /** Convergence tolerance */
  private static final double TOLERANCE = 1e-2;
  
  /** Number of latent factors */
  protected static final int NLATENT = 20;
  
  private AlsGraph mGraph;
  
  protected AlsUpdater(AlsGraph graph) {
    mGraph = graph;
  }

  @Override
  protected AlsUpdater clone() {
    return new AlsUpdater(mGraph);
  }

  @Override
  public void update(Context context, AlsVertex vertex) {
  
    vertex.mSSE = 0;
  
    // if there are no neighbors just return -----------------------------------
    Set<DefaultWeightedEdge> edges = mGraph.edgesOf(vertex);
    if (edges.isEmpty()) return;
  
    // compute X ---------------------------------------------------------------
    DenseMatrix X = new DenseMatrix(edges.size(), NLATENT);
    DenseVector y = new DenseVector(edges.size());
    
    int i=0;
    for (final DefaultWeightedEdge edge : edges) {
      // set x values
      Vector neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex).vector();
      for (int j = 0; j < NLATENT; j++) X.set(i, j, neighbor.get(j));
      // set rating
      y.set(i, mGraph.getEdgeWeight(edge));
      i++;
    }
    
    // compute X'X and X'y -----------------------------------------------------
    DenseMatrix XtX = new DenseMatrix(NLATENT, NLATENT);
    X.transAmult(X, XtX);
    DenseVector Xty = new DenseVector(NLATENT);
    X.transMult(y, Xty); 
  
    // regularization
    for (i = 0; i < NLATENT; i++) XtX.add(i, i, (LAMBDA) * edges.size());
  
    // solve the least squares problem -----------------------------------------
    double[] weights = DenseCholesky.factorize(XtX).solve(new DenseMatrix(Xty)).getData();
    vertex.setVector(new DenseVector(weights));
  
    // update the RMSE and reschedule neighbors --------------------------------
    for (final DefaultWeightedEdge edge : edges) {
  
      // get the neighbor id
      AlsVertex neighbor = Graphs.getOppositeVertex(mGraph, edge, vertex);
      final double pred = vertex.vector().dot(neighbor.vector());
      final double error = Math.abs(mGraph.getEdgeWeight(edge) - pred);
      vertex.mSSE += error * error;
  
      // reschedule neighbors
      if (error > TOLERANCE && vertex.mResidual > TOLERANCE)
        context.schedule(neighbor, new AlsUpdater(mGraph));
  
    }
  
  } // end of operator()

}
