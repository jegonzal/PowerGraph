package org.graphlab.data;

public class ScalarVertex implements Vertex {

  /** GraphLab (or proxy vertex) ID */
  private int mRawId = 0;
  
  /** Application ID */
  private int mId;
  
  private double mValue;

  /**
   * @param id
   *      ID used in TSV-formatted files
   */
  public ScalarVertex(int id){
    mId = id;
  }
  
  public int id(){
    return mId;
  }
  
  public int rawId() {
    return mRawId;
  }

  public void setRawId(int id) {
    mRawId = id;
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
  
  @Override
  public int hashCode(){
    return mId;
  }
  
  @Override
  public boolean equals(Object other){
    if (!(other instanceof ScalarVertex)) return false;
    ScalarVertex otherVertex = (ScalarVertex) other;
    return mId == otherVertex.mId;
  }

}
