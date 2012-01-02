package org.graphlab;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.graphlab.data.Edge;
import org.graphlab.data.Graph;
import org.graphlab.data.Vertex;

/**
 * GraphLab Core. This interfaces with the C++ library via JNI.
 * 
 * @author Jiunn Haur Lim
 * @param <G> 		graph type -> must extend {@link org.graphlab.data.Graph}
 */
public final class Core<G extends Graph<? extends Vertex, ? extends Edge>> {

	/** Logger (Java layer only, the C++ layer has its own logger.) */
	private static final Logger logger = Logger.getLogger (Core.class);

	/** Indicates if the core has been destroyed and cannot be reused. */
	private boolean mDestroyed = false;

	/** Address of graphlab::core */
	private long mCorePtr;

	/** Mapping from application vertex IDs to graphlab vertex IDs */
	private Map<Integer, Integer> mIdMap = null;
	
	/** Updaters scheduled by {@link #schedule(int, Updater)} */
	private List<Updater> mUpdaters;

	static {
		// load the JNI library
		System.loadLibrary("graphlabjni");
		logger.trace ("JNI library libgraphlabjni loaded.");
	}

	/**
	 * Creates a new GraphLab core.
	 * <b>Call {@link #destroy()} when done to free up resources.</b>
	 * 
	 * @throws CoreException if there was an error creating the core
	 */
	public Core() throws CoreException {

		mCorePtr = createCore();

		if (0 >= mCorePtr)
			throw new CoreException("Unable to create a core.");

		logger.trace ("Core created.");
		
		mUpdaters = Collections.synchronizedList(new ArrayList<Updater>());

	}

	/**
	 * Tells core to operate on this graph. This creates a proxy graph in the
	 * GraphLab engine.
	 * 
	 * Once this is called, your graph must not be modified (or there be dragons!).
	 * 
	 * @param graph
	 *            the graph to operate on
	 */
	public void setGraph(G graph) {

		if (null == graph)
			throw new NullPointerException("graph must not be null.");
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");

		// inspect graph size
		int size = graph.size();
		if (size <= 0)
			throw new IllegalArgumentException("graph must not be empty ");
		
		// inspect vertices
		Collection<? extends Vertex> vertices = graph.vertices();
		if (null == vertices)
			throw new IllegalArgumentException("graph must not be empty.");
		
		long startTime = System.currentTimeMillis();
		
		// create new map from application vertex IDs to graphlab vertex IDs
		mIdMap = new HashMap<Integer, Integer>(size);

		// add vertices
		for (Vertex v : vertices) {
			mIdMap.put(v.id(), addVertex(mCorePtr, v.id()));
		}
		Context.getInstance().setIdMap(mIdMap);

		// add edges
		for (Vertex v : vertices) {
			for (Edge e : graph.outgoingEdges(v.id())) {
				addEdge(mCorePtr, mIdMap.get(e.source()),
						mIdMap.get(e.target()));
			}
		}

		long elapsed = System.currentTimeMillis() - startTime;
		logger.info ("Graph transferring took: " + elapsed + " ms.");

	}

	/**
	 * Destroys the GraphLab core. Once destroyed, this object may not be used.
	 */
	public void destroy() {

		if (mDestroyed) {
			logger.warn("Core has already been destroyed and may not be destroyed again.");
			return;
		}

		if (0 == mCorePtr)
			throw new IllegalStateException("Core missing or was never allocated.");

		destroyCore(mCorePtr);
		mDestroyed = true;
		
		// remove references to allow garbage collection
		mIdMap = null;
		mUpdaters = null;

		logger.trace("Core destroyed.");

	}

	/**
	 * Run the engine until a termination condition is reached or there are no
	 * more tasks remaining to execute.
	 * 
	 * @return run time
	 */
	public double start() {
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");
		
		logger.info("GraphLab engine started.");
		return start(mCorePtr);
		
	}

	/**
	 * Schedule the execution of an update functor on a particular vertex.
	 * 
	 * @param vertexId
	 *            application vertex ID
	 * @param updater
	 * 			  updater to execute
	 * @throws NoSuchElementException
	 *             if the specified vertex did not exist in the graph that was
	 *             passed to {@link #setGraph(Graph)}.
	 */
	public void schedule(int vertexId, Updater updater) {

		if (null == updater)
			throw new NullPointerException("updater must not be null.");
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");

		Integer id = mIdMap.get(vertexId);
		if (null == id)
			throw new NoSuchElementException("vertex did not exist in the graph that was passed to #setGraph.");

		if (updater.id() != Updater.ID_NOT_SET){
			// this updater already has an ID, insert
			mUpdaters.set(updater.id(), updater);
		}else {
			// otherwise, add to end of list and assign id
			mUpdaters.add(updater);
			updater.setId(mUpdaters.size()-1);
		}
		
		schedule(mCorePtr, id, updater.id());
		
	}

	/**
	 * Executes the updater on the specified vertex. This is <em>only</em>
	 * invoked by the proxy updater in the JNI library.
	 * 
	 * @param contextPtr
	 *         address of C++ context object
	 * @param vertexId
	 * 				application vertex ID
	 * @param updaterId
	 * 				updater ID (as assigned by {@link #schedule(int, Updater)}).
	 */
	private void execUpdate (long corePtr, long contextPtr, int vertexId, int updaterId){
		
		Updater updater = mUpdaters.get(updaterId);
		updater.update(corePtr, contextPtr, vertexId);
		
	}
	
	/**
	 * Creates and initializes graphlab::core -> dynamically allocates a core.
	 * Must be freed by a corresponding call to {@link #destroyCore()}.
	 * 
	 * @return address of core or 0 on failure
	 */
	private native long createCore();

	/**
	 * Deletes the graphlab::core that was allocated in {@link #initCore()}.
	 * 
	 * @param ptr
	 *            address of core
	 */
	private native void destroyCore(long ptr);
	
	/**
	 * Add additional vertices up to <tt>n</tt>. This will fail if resizing down. 
	 * @param ptr
	 * 			{@link #mCorePtr}
	 * @param n
	 * 			number of vertices to add
	 */
	private native void resizeGraph(long ptr, int n);
	
	/**
	 * Adds a vertex to the native graph.
	 * @param ptr
	 * 			{@link #mCorePtr}
	 * @param id
	 * 			application vertex ID
	 * @return
	 * 			graphlab vertex ID
	 */
	private native int addVertex(long ptr, int id);

	/**
	 * Adds an edge to the native graph.
	 * @param ptr
	 * 			{@link #mCorePtr}
	 * @param source
	 * 			graphlab vertex ID
	 * @param target
	 * 			graphlab vertex ID
	 */
	private native void addEdge(long ptr, int source, int target);
	
	/**
	 * Add a single update function to a single vertex
	 * @param ptr
	 * @param vertex_id
	 * @param updater_id
	 */
	private native void schedule(long ptr, int vertex_id, int updater_id);

	/**
	 * Run the engine until a termination condition is reached or there are no more tasks remaining to execute. 
	 * @param ptr
	 * @return runtime
	 */
	private native double start(long ptr);

	/**
	 * Generic exception for dealing with GraphLab cores
	 * 
	 * @author Jiunn Haur Lim
	 */
	public static class CoreException extends Exception {

		private static final long serialVersionUID = -1231034883125084972L;

		public CoreException(String string) {
			super(string);
		}

	}

}
