package org.graphlab.data;

public interface Edge {
	
	/**
	 * @return ID of source vertex
	 */
	public int source();
	
	/**
	 * @return of target vertex
	 */
	public int target();

}
