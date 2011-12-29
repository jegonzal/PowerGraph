package org.graphlab.data;

/**
 * Implement this interface to store more data in
 * the edges.
 * @author Jiunn Haur Lim
 */
public interface Edge {
	
	/**
	 * @return ID of source vertex
	 */
	public int source();
	
	/**
	 * @return ID of target vertex
	 */
	public int target();

}
