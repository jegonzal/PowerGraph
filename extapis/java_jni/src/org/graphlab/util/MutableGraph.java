package org.graphlab.util;

/**
 * For {@link org.graphlab.util.GraphLoader} only. The GraphLab engine
 * doesn't need this.
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public interface MutableGraph<V,E> {

	public int addVertex (V vertex);
	public void addEdge (E edge);
	public int size();
	
}
