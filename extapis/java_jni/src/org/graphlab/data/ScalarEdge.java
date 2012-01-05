package org.graphlab.data;

/**
 * Simple weighted edge.
 * 
 * @author akyrola
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ScalarEdge implements Edge {

  private int source;
  private int target;
  private double weight;

  public ScalarEdge(int source, int target, double weight) {
    this.source = source;
    this.target = target;
    this.weight = weight;
  }

  public int source() {
    return source;
  }

  public int target() {
    return target;
  }

  public double weight() {
    return weight;
  }

  public void setWeight(double weight) {
    this.weight = weight;
  }

}
