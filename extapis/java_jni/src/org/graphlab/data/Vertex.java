package org.graphlab.data;


/**
 * Vertex for graphs in GraphLab. All graphs that are passed to
 * {@link org.graphlab.Core} must use this vertex type.
 * 
 * <p>
 * When {@link org.graphlab.Core#setGraph(org.jgrapht.DirectedGraph)} is called,
 * a proxy C++ graph is created, and a proxy vertex is created for each Vertex
 * object. {@link #setRawId(int)} is called to save the ID of the proxy vertex
 * in the corresponding Vertex object.
 * {@link org.graphlab.Context#schedule(Vertex, org.graphlab.Updater)} and
 * {@link org.graphlab.Core#schedule(Vertex, org.graphlab.Updater)} needs to
 * know the ID of the proxy vertex in order to schedule updates.
 * </p>
 * 
 * <p>
 * Because we are using JGraphT, your Vertex implementations must also override
 * {@link #hashCode()} and {@link #equals(Object)}.
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @see <a href="http://www.jgrapht.org/">JGraphT</a>
 */
public interface Vertex {

	/**
	 * @return ID of corresponding proxy vertex.
	 */
	public int rawId();

	/**
	 * Sets ID of corresponding proxy vertex. This is done during
	 * {@link org.graphlab.Core#setGraph(org.jgrapht.DirectedGraph)}.
	 * Applications should <em>not</em> invoke this method.
	 * @param id
	 *            id of this vertex's proxy vertex.
	 */
	public void setRawId(int id);

}
