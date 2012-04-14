package org.graphlab;

/**
 * Collection of built-in consistency models.
 * 
 * <p>The scope of a vertex is the data on the vertex, the data on all adjacent
 * edges as well as data on all adjacent vertices. An update function executed
 * on a vertex has the ability to read and write to all data within the scope of
 * the vertex.</p>
 * 
 * <p>The ability for an update function to read or write a piece of data does not
 * mean that it can do so in a consistent fashion. For instance, without
 * adequate protection, it is conceivable that update functions executed on
 * neighboring vertices could modify the same piece of data at the same time.</p>
 * 
 * <p>GraphLab therefore provides the concept of a consistency class which
 * expresses both the level of protection provided, as well as the amount of
 * parallelism available. The consistency type can be selected through
 * {@link CoreConfiguration#setScopeType(Scope)} and {@link Core#setScopeType(Scope)}.</p>
 * 
 * <p>There are three basic types of scopes (as well a number of undocumented ones
 * which you probably should not use). The scopes are currently implemented
 * through the use of a read-write lock on each vertex. Different types of
 * scopes then make use of different locking schemes.</p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://graphlab.org/doxygen/html/Scopes.html">Scopes</a>
 */
public enum Scope {

  /**
   * This scope type is the safest and guarantees safe access to the current
   * vertex, adjacent edges, and adjacent vertices. This is implemented by
   * acquiring a write lock on all vertices in the scope.
   */
  FULL("full"),

  /**
   * This scope type guarantees safe read/write access to the current vertex as
   * well as adjacent edges. In addition, you can also read (but not write) data
   * on adjacent vertices safely. This is implemented by acquiring a write lock
   * on the current vertex, and a read lock on all adjacent vertices.
   */
  EDGE("edge"),

  /**
   * This scope type guarantees only safe read/write access to the current
   * vertex. This is implemented by acquiring a write lock on the current
   * vertex.
   */
  VERTEX("vertex");

  private final String mStr;

  private Scope(String str) {
    mStr = str;
  }

  @Override
  public String toString() {
    return mStr;
  }

}
