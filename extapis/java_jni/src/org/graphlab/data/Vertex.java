package org.graphlab.data;

/**
 * Vertex for graphs in GraphLab. All graphs that are passed to
 * {@link org.graphlab.Core} must use this vertex type.
 * 
 * <p>
 * When {@link org.graphlab.Core#setGraph(org.jgrapht.DirectedGraph)} is called,
 * a proxy C++ graph is created and a proxy vertex is created for each Vertex
 * object. {@link #setRawId(int)} is called to save the ID of the proxy vertex
 * in the corresponding Vertex object.
 * {@link org.graphlab.Context#schedule(Vertex, org.graphlab.Updater)} and
 * {@link org.graphlab.Core#schedule(Vertex, org.graphlab.Updater)} needs to
 * know the ID of the proxy vertex in order to schedule updates.
 * </p>
 * 
 * <p>
 * On top of rawId, each vertex also needs an id (we call this the application vertex
 * id) so that JGraphT can uniquely identify the vertex. The default implementation
 * assumes that two vertices are equal if and only if they share the same application
 * vertex id. Override {@link #hashCode()} and {@link #equals(Object)} if you need
 * more sophisticated behavior.
 * </p>
 * 
 * <p>Note to self: this is all too clunky. After Joey makes the changes to the
 * core graphs, merge application vertex ID and raw ID.</p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://www.jgrapht.org/">JGraphT</a>
 */
public abstract class Vertex {

  /** GraphLab (or proxy vertex) ID */
  private int mRawId = 0;
  
	/**
	 * @return ID of corresponding proxy vertex.
	 */
  public final int rawId() {
    return mRawId;
  }

	/**
	 * Sets ID of corresponding proxy vertex. This is done during
	 * {@link org.graphlab.Core#setGraph(org.jgrapht.DirectedGraph)}.
	 * Applications should <em>not</em> invoke this method.
	 * @param id
	 *            id of this vertex's proxy vertex.
	 */
  public final void setRawId(int id) {
    mRawId = id;
  }
  
  /**
   * Specifies a vertex ID. Required by JGraphT
   * to differentiate between vertices.
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
