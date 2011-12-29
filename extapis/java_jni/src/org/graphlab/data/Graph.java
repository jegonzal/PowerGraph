package org.graphlab.data;

import java.util.Collection;

/**
 * Graph interface.
 * 
 * This is the interface between GraphLab and your application graphs.
 * 
 * GraphLab graphs are directional. Self-edges and duplicate edges are not
 * allowed. In this new graph interface, only the graph object is responsible
 * for maintaining structure; the vertices and edges are only responsible for
 * maintaining data (although edge will know its source and target.) This is
 * done to mirror graphlab::graph more closely (so that future changes would
 * hopefully be less drastic), reduce data duplication, and allow more graph
 * implementations e.g. adjacency matrices.
 * 
 * Application graphs must implement this interface in order to be passed to the
 * {@link org.graphlab.Core}. This interface was designed to be minimal so that
 * it could be pasted easily over any existing graphs.
 * 
 * @param <V>
 *            vertex type
 * @param <E>
 *            edge type
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public interface Graph<V extends Vertex, E extends Edge> {

	/**
	 * Obtains an immutable collection of all the vertices in this graph.
	 * 
	 * Implementations should ensure that this collection is has no duplicate
	 * vertices (for performance purposes, this method only enforces a
	 * Collection interface and not a Set interface - which means that you could
	 * just use an array).
	 * 
	 * @return collection of vertices
	 */
	public Collection<V> vertices();

	/**
	 * Retrieves a vertex given its ID.
	 * 
	 * @param id
	 *            id of vertex (see {@link Vertex#id()}).
	 * @return corresponding vertex if found, or null
	 */
	public V getVertex(int id);

	/**
	 * Retrieves an immutable collection of outgoing edges from this vertex.'
	 * 
	 * Implementations should ensure that this collection is has no duplicate
	 * edges (for performance purposes, this method only enforces a Collection
	 * interface and not a Set interface - which means that you could just use
	 * an array).
	 * 
	 * @param id
	 * @return outgoing edges
	 * @throws NoSuchElementException
	 * 				if the vertex is not found
	 */
	public Collection<E> outgoingEdges(int id);

	/**
	 * @return number of vertices in the graph
	 */
	public int size();

}
