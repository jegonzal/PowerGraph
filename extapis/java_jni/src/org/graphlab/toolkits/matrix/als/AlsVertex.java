package org.graphlab.toolkits.matrix.als;

import org.graphlab.toolkits.matrix.VectorVertex;

/**
 * Vertex for ALS.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class AlsVertex extends VectorVertex {

  public double mResidual = Double.MAX_VALUE;
  public int mNUpdates = 0;
  public double mSquaredError = Double.MAX_VALUE;
  
  public AlsVertex(){
    super();
  }
  
  public AlsVertex(int id) {
    super(id);
  }

}
