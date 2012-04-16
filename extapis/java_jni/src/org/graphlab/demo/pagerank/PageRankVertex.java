package org.graphlab.demo.pagerank;

import org.graphlab.data.ScalarVertex;

/**
 * PageRankVertex.
 * 
 * Contains some additional data.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class PageRankVertex extends ScalarVertex {

  public int mNUpdates = 0;
  public double mOldValue = 0.0;
  
  public PageRankVertex (){
    // set default value to 1
    setValue(1);
  }
  
}
