package org.graphlab.data;

/**
 * Vertex that holds a single value
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ScalarVertex extends Vertex {
  
  /** Application ID */
  private int mId;
  
  private double mValue;

  /**
   * The default constructor does nothing. Make sure you invoke
   * {@link ScalarVertex#setId}.
   */
  public ScalarVertex(){}
  
  /**
   * @param id
   *      ID used in TSV-formatted files. Required by JGraphT
   *      to differentiate between vertices.
   */
  public ScalarVertex(int id){
    this(id, 0);
  }
  
  /**
   * @param id
   *      ID used in TSV-formatted files. Required by JGraphT
   *      to differentiate between vertices.
   * @param value
   *      initial value of vertex
   */
  public ScalarVertex(int id, double value){
    mId = id;
    mValue = value;
  }
  
  public int id(){
    return mId;
  }
  
  /**
   * Specifies a vertex ID. Required by JGraphT
   * to differentiate between vertices.
   * @param id
   */
  public void setId(int id){
    mId = id;
  }
  
  /**
   * @return the mValue
   */
  public double value() {
    return mValue;
  }

  /**
   * @param d the mValue to set
   */
  public void setValue(double d) {
    this.mValue = d;
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
    if (!(other instanceof ScalarVertex)) return false;
    ScalarVertex otherVertex = (ScalarVertex) other;
    return mId == otherVertex.mId;
  }
  
  @Override
  public String toString(){
    return "id: " + mId + " value: " + mValue;
  }

}
