package org.graphlab.toolkits.matrix;

import org.graphlab.data.Vertex;

import cern.colt.matrix.AbstractMatrix1D;

/**
 * Vertex that holds a vector.
 * 
 * In the context of a recommender system, a vertex is either a user
 * or an object of interest.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class VectorVertex extends Vertex {

  /**
   * Application vertex ID.
   * In the context of a recommender system, this will be the ID
   * of the user and the ID of the object.
   */
  private int mId;
  
  /** Data for this vertex */
  private AbstractMatrix1D mVector;
  
  /**
   * Constructs a new vector vertex with the given ID
   * @param id
   *          Used to identify this vertex within the application.
   */
  public VectorVertex(int id){
    mId = id;
  }
  
  public int id(){
    return mId;
  }
  
  /**
   * Sets the data for this vertex.
   * @param vector
   */
  public void setVector(AbstractMatrix1D vector){
    mVector = vector;
  }
  
  /*
   * IMPORTANT: must override this for JGraphT
   * (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode(){
    return mId;
  }
  
  /*
   * IMPORTANT: must override this for JGraphT
   * (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object other){
    if (!(other instanceof VectorVertex)) return false;
    return mId == ((VectorVertex) other).mId;
  }
  
  @Override
  public String toString(){
    return "id: " + mId + " value: " + mVector;
  }
  
}
