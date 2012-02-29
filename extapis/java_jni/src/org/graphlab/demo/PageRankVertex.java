package org.graphlab.demo;

import org.graphlab.data.ScalarVertex;

/**
 * PageRankVertex.
 * 
 * Contains some additional data.
 * @author Jiunn Haur Lim
 *
 */
public class PageRankVertex extends ScalarVertex {

  public int mNUpdates = 0;
  public double mOldValue = 0.0;
  
}
