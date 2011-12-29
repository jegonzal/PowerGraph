package org.graphlab;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.logging.Logger;

import org.graphlab.data.Edge;
import org.graphlab.data.Graph;
import org.graphlab.data.Vertex;

/**
 * GraphLab Core. This interfaces with the C++ library via JNI.
 * 
 * @author Jiunn Haur Lim
 * @param <G> 		graph type -> must extend {@link org.graphlab.data.Graph}
 */
public class Core<G extends Graph<? extends Vertex, ? extends Edge>> {

	/** Name of logger that this class uses */
	public static final String TAG = "org.graphlab.core";

	/** Logger (Java layer only, the C++ layer has its own logger.) */
	private static final Logger logger = Logger.getLogger(TAG);

	/** Indicates if the core has been destroyed */
	private boolean mDestroyed = false;

	/** Address of graphlab::core */
	private long mCorePtr;

	/** Mapping from application vertex IDs to graphlab vertex IDs */
	private Map<Integer, Integer> mIdMap = null;
	
	private List<Updater> mUpdaters;

	static {
		// load the JNI library
		System.loadLibrary("graphlabjni");
		logger.fine("graphlabjni loaded.");
	}

	/**
	 * Creates a new GraphLab core. Call {@link #destroy()} when done to free up
	 * resources.
	 * 
	 * @throws Exception
	 */
	public Core() throws Exception {

		mCorePtr = createCore();

		if (0 == mCorePtr) {
			// TODO
			throw new CoreException("Unable to create a core.");
		}

		logger.finer("Core created.");
		mUpdaters = Collections.synchronizedList(new ArrayList<Updater>());

	}

	/**
	 * Tells core to operate on this graph. Once this is called, your graph must
	 * not be modified (or there be dragons!)
	 * 
	 * @param graph
	 *            the graph to operate on
	 */
	public void setGraph(G graph) {

		if (null == graph)
			throw new NullPointerException("graph must not be null.");

		long startTime = System.currentTimeMillis();

		int size = graph.size();

		mIdMap = new HashMap<Integer, Integer>(size);

		// add vertices
		for (Vertex v : graph.vertices()) {
			mIdMap.put(v.id(), addVertex(mCorePtr, v.id()));
		}

		// add edges
		for (Vertex v : graph.vertices()) {
			for (Edge e : graph.outgoingEdges(v.id())) {
				addEdge(mCorePtr, mIdMap.get(e.source()),
						mIdMap.get(e.target()));
			}
		}

		long elapsed = System.currentTimeMillis() - startTime;
		logger.fine("Graph transferring took: " + elapsed + " ms.");

	}

	/**
	 * Destroys the GraphLab core. Once destroyed, this object may not be used.
	 */
	public void destroy() {

		if (mDestroyed) {
			logger.warning("Core has already been destroyed and may not be destroyed again.");
			return;
		}

		if (0 == mCorePtr)
			throw new IllegalStateException(
					"Core missing or was never allocated.");

		destroyCore(mCorePtr);
		mDestroyed = true;
		mIdMap = null;

		logger.finer("Core destroyed.");

	}

	/**
	 * Run the engine until a termination condition is reached or there are no
	 * more tasks remaining to execute.
	 * 
	 * @return run time
	 */
	public double start() {
		return start(mCorePtr);
	}

	/**
	 * Schedule the execution of an update functor on a particular vertex.
	 * 
	 * @param vertex_id
	 *            application vertex ID
	 * @param updater
	 * @throws NoSuchElementException
	 *             if the specified vertex did not exist in the graph that was
	 *             passed to {@link #setGraph(Graph)}.
	 */
	public void schedule(int vertex_id, Updater updater) {

		if (null == updater)
			throw new NullPointerException("updater must not be null.");

		Integer id = mIdMap.get(vertex_id);
		if (null == mIdMap.get(vertex_id))
			throw new NoSuchElementException(
					"vertex did not exist in the graph that was passed to #setGraph.");

		if (updater.id() != Updater.ID_NOT_SET){
			// this updater already has an ID
			mUpdaters.set(updater.id(), updater);
		}else {
			mUpdaters.add(updater);
			updater.setId(mUpdaters.size()-1);
		}
		
		schedule(mCorePtr, id, updater.id());
		
	}

	public void execUpdate (int vertex_id, int updater_id){
		
		Updater updater = mUpdaters.get(updater_id);
		updater.update(vertex_id);
		
	}
	
	// TODO: documentation and error handling
	private native void schedule(long ptr, int vertex_id, int updater_id);

	/**
	 * Add additional vertices up to provided num_vertices. This will fail if resizing down. 
	 * @param ptr
	 * 			mCorePtr
	 * @param n
	 * 			number of vertices to add.
	 * TODO: error handling
	 */
	private native void resizeGraph(long ptr, int n);

	// TODO: documentation and error handling
	private native int addVertex(long ptr, int id);

	// TODO: documentation and error handling
	private native void addEdge(long ptr, int source, int target);

	// TODO: documentation and error handling
	private native double start(long ptr);

	/**
	 * Creates and initializes graphlab::core -> dynamically allocates a core.
	 * Must be freed by a corresponding call to {@link #destroyCore()}.
	 * 
	 * @return address of core
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
