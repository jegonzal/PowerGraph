package org.graphlab.data;

/**
 * Vertex for graphs in GraphLab. Because a proxy C++ graph is maintained,
 * every vertex needs to know the ID of its corresponding proxy vertex.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public interface Vertex {

	/**
	 * @return ID of corresponding proxy vertex.
	 */
	public int rawId();

	/**
	 * Sets ID of corresponding proxy vertex. This is done during
	 * {@link org.graphlab.Core#setGraph(org.jgrapht.DirectedGraph)}.
	 * @param id
	 *            id of this vertex's proxy vertex.
	 */
	public void setRawId(int id);

}
