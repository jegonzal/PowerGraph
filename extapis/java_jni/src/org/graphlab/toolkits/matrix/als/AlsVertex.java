package org.graphlab.toolkits.matrix.als;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.graphlab.toolkits.matrix.VectorVertex;

/**
 * Vertex for ALS. Contains residual, squared error, latent factors, and 
 * biases.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class AlsVertex extends VectorVertex {
  
  /**
   * Application vertex ID.
   * In the context of a recommender system, this will be the ID
   * of the user and the ID of the object.
   */
  private int mId;

  /** Average change in latent factors */
  public double mResidual = Double.POSITIVE_INFINITY;
  
  /** Sum of squared errors */
  public double mSSE = Double.POSITIVE_INFINITY;
  
  public AlsVertex(){
    super();
  }
  
  public AlsVertex(int id) {
    mId = id;
  }
  
  @Override
  public int id(){
    return mId;
  }
  
  @Override
  public void setId(int id) {
    mId = id;    
  }
  
  /**
   * Update vectors and calculate average change in latent factors
   */
  @Override
  public void setVector(Vector vector){
    Vector oldVector = vector();
    super.setVector(vector);
    if (null == oldVector) oldVector = new DenseVector(vector.size());
    Vector difference = oldVector.add(-1.0, vector); // careful - this overwrites oldVector
    mResidual = difference.norm(Vector.Norm.One) / difference.size();
  }

  @Override
  public String toString(){
    return "id: " + mId + " value: " + vector();
  }

}
