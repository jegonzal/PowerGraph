package org.graphlab.data;

/**
 * Implement this interface to store more data in the vertices.
 * @author Jiunn Haur Lim
 */
public interface Vertex {

	/**
	 * ID of this vertex. Must uniquely identify the vertex within a graph. The
	 * GraphLab core will use this ID to forward update calls from the graph
	 * maintained by the engine to your application graph.
	 * 
	 * @return id of this vertex
	 */
	// TODO: other data type for id?
	public int id();

	/**
	 * Sets ID of this vertex. Must uniquely identify the vertex within a graph.
	 * The GraphLab core will use this ID to forward update calls from the graph
	 * maintained by the engine to your application graph.
	 * 
	 * @param id
	 *            id of this vertex
	 */
	public void setId(int id);

}
