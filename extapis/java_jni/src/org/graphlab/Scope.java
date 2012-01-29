package org.graphlab;

/**
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://graphlab.org/doxygen/html/Scopes.html">Scopes</a>
 */
public enum Scope {

  /** This ensures full data consistency within the scope */
  FULL("full"),

  /** This ensures data consistency with just the vertex and edges */
  EDGE("edge"),

  /**
   * This ensures that a vertex cannot be updated by two processors
   * simultaneously
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
