package org.graphlab.data;

/**
 * Vertex for graphs in GraphLab. All graphs that are passed to
 * {@link org.graphlab.Core} must use this vertex type.
 * 
 * <p>
 * When {@link org.graphlab.Core#setGraph} is called, a proxy C++ graph is
 * created and a proxy vertex is created for each Vertex object. {@link #id()}
 * is used to uniquely identify/reference vertices.
 * </p>
 * 
 * <p>
 * The default implementation assumes that two vertices are equal if and only if
 * they share the same vertex id. Override {@link #hashCode()} and
 * {@link #equals(Object)} if you need more sophisticated behavior.
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://www.jgrapht.org/">JGraphT</a>
 */
public abstract class Vertex {
  
  /**
   * Specifies a vertex ID. Required by JGraphT to differentiate between
   * vertices.
   * 
   * @param id
   */
  public abstract void setId(int id);
  
  public abstract int id();
  
  /*
   * IMPORTANT: must override this for JGraphT
   * (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode(){
    return id();
  }
  
  /*
   * IMPORTANT: must override this for JGraphT
   * (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object other){
    if (this == other) return true;
    if (!(other instanceof Vertex)) return false;
    return id() == ((Vertex) other).id();
  }

}
