package org.graphlab.data;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

import org.graphlab.util.MutableGraph;

/**
 * Efficient representation for sparse graphs. Default Graph class to be used.
 * 
 * @author akyrola
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class SparseGraph<V extends Vertex, E extends Edge>
	implements Graph<V, E>, MutableGraph<V, E> {

	/** Collection of vertices */
	private V[] mVertices;

	/**
	 * Collection of outgoing edges. Since it's sparse, each vertex is likely to have a
	 * small number of edges.
	 */
	private List<E>[] mOutgoingEdges;
	
	/**
	 * Collection of incoming edges. Since it's sparse, each vertex is likely to have a
	 * small number of edges.
	 */
	private List<E>[] mIncomingEdges;

	/** Number of vertices */
	private int mSize;

	/** Default estimated number of vertices */
	public static final int DEFAULT_ESTIMATE = 1000;

	/**
	 * Constructor
	 */
	public SparseGraph() {
		this(DEFAULT_ESTIMATE);
	}

	/**
	 * Creates a new graph with estimated number of vertices. Internal data
	 * structures will grow if the estimate is exceeded.
	 * 
	 * @param n
	 *            estimated number of vertices
	 */
	@SuppressWarnings("unchecked")
	public SparseGraph(int n) {
		mVertices = (V[]) new Vertex[n];
		mOutgoingEdges = (List<E>[]) new List[n];
		mIncomingEdges = (List<E>[]) new List[n];
		mSize = 0;
	}

	@Override
	public V getVertex(int id) {
		if (!isVertexInGraph(id))
			return null;
		return mVertices[id];
	}

	/**
	 * Adds a vertex to the graph and assigns an ID to it.
	 * 
	 * @param vertex
	 *            the vertex to be added
	 * @return the assigned vertex id
	 */
	@Override
	public synchronized int addVertex(V vertex) {
		if (null == vertex)
			throw new NullPointerException("vertex must not be null.");
		ensureCapacity(mSize + 1);
		mVertices[mSize++] = vertex;
		vertex.setId(mSize - 1);
		return vertex.id();
	}

	// TODO: sort edges and use binary search
	/**
	 * Returns an edge from source vertex to target vertex. Note: this is a slow
	 * operation (O(degree)) since edges are not sorted.
	 * 
	 * @param source
	 *            vertex id of the origin vertex
	 * @param target
	 *            vertex id of the destination vertex
	 * @return edge spanning source to target vertex
	 */
	public Edge getEdge(int source, int target) {
		if (!isVertexInGraph(source) || !isVertexInGraph(target))
			return null;
		for (E edge : mOutgoingEdges[source])
			if (edge.target() == target)
				return edge;
		return null;
	}

	/**
	 * Adds an edge between two vertices in the graph
	 * 
	 * @param edge
	 *            the edge to be added
	 */
	@Override
	public synchronized void addEdge(E edge) {

		if (null == edge)
			throw new NullPointerException("e must not be null.");

		// ensure that source exists
		int source = edge.source();
		if (!isVertexInGraph(source))
			throw new NoSuchElementException("source not found in graph.");

		// ensure that target exists
		int target = edge.target();
		if (!isVertexInGraph(target))
			throw new NoSuchElementException("target not found in graph.");

		if (null == mOutgoingEdges[source])
			mOutgoingEdges[source] = new LinkedList<E>();
		mOutgoingEdges[source].add(edge);
		
		if (null == mIncomingEdges[target])
			mIncomingEdges[target] = new LinkedList<E>();
		mIncomingEdges[target].add(edge);

	}

	@Override
	public Collection<V> vertices() {
		// TODO: optimize and remove copy (instead, use a proxy)
		return Collections.unmodifiableList(
			Arrays.asList(Arrays.copyOf(mVertices, mSize))
		);
	}

	@Override
	public int size() {
		return mSize;
	}

	@Override
	public Collection<E> outgoingEdges (int id){
		if (!isVertexInGraph(id))
			throw new NoSuchElementException ("vertex " + id + " not found in this graph.");
		if (null == mOutgoingEdges[id]) return new LinkedList<E>();
		return Collections.unmodifiableList(mOutgoingEdges[id]);
	}
	
	public Collection<E> incomingEdges (int id){
		if (!isVertexInGraph(id))
			throw new NoSuchElementException ("vertex " + id + " not found in this graph.");
		if (null == mIncomingEdges[id]) return new LinkedList<E>();
		return Collections.unmodifiableList(mIncomingEdges[id]);
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer(size() * 30);
		for (Vertex v : mVertices) {
			if (v != null) {
				sb.append(v.id());
				sb.append(", ");
				sb.append(v);
				sb.append(", ");
				sb.append(outgoingEdges(v.id()));
				sb.append("\n");
			}
		}
		return sb.toString();
	}
	
	/**
	 * Checks if the vertex corresponding to the given ID exists in this graph.
	 * 
	 * @param id
	 * @return true if exists, false otherwise
	 */
	private boolean isVertexInGraph(int id) {
		return !(id < 0 || id >= mSize);
	}

	private void ensureCapacity(int size) {

		// if large enough, just return
		if (mVertices.length >= size)
			return;

		// grow by 100% for constant amortized cost
		mVertices = Arrays.copyOf(mVertices, mVertices.length * 2);
		mOutgoingEdges = Arrays.copyOf(mOutgoingEdges, mOutgoingEdges.length * 2);
		mIncomingEdges = Arrays.copyOf(mIncomingEdges, mIncomingEdges.length * 2);

	}

}
