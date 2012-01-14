package org.graphlab;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.graphlab.data.Edge;
import org.graphlab.data.Graph;
import org.graphlab.data.Vertex;

/**
 * GraphLab Core.
 * 
 * <p>
 * This interfaces with the C++ library via
 * <abbr title="Java Native Interface">JNI</abbr> and mirrors
 * <tt>graphlab::core</tt>.
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 * @param <G>
 *          graph type; must extend {@link org.graphlab.data.Graph}
 */
public final class Core<G extends Graph<? extends Vertex, ? extends Edge>> {

	/** Logger (Java layer only, the C++ layer has its own logger) */
	private static final Logger logger = Logger.getLogger (Core.class);

	/** Indicates if the core has been destroyed and cannot be reused. */
	private boolean mDestroyed = false;

	/** Address of graphlab::core object */
	private long mCorePtr;
	
  /** Mapping from application vertex IDs to graphlab vertex IDs */
  private Map<Integer, Integer> mIdMap = null;

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

	}

	/**
	 * Tells core to operate on this graph. This creates a proxy graph in the
	 * GraphLab engine.
	 * 
	 * <p>
	 * Once this is called, your graph must not be modified (or there be dragons!)
	 * </p>
	 * 
	 * @param graph
	 *            the graph to operate on
	 * @throws NullPointerException if graph is null
	 * @throws IllegalArgumentException if graph is empty
	 * @throws IllegalStateException if {@link #destroy()} was invoked on this object
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
   * Calling {@link #destroy()} more than once on the same object will generate
   * a warning in the logs but will not have an effect.
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
		logger.trace("Core destroyed.");
		
		// remove references to allow garbage collection
		mIdMap = null;

	}

  /**
   * Run the engine until a termination condition is reached or there are no
   * more tasks remaining to execute.
   * 
   * @return run time
   * @throws IllegalStateException if {@link #destroy()} was invoked on this object
   */
	public double start() {
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");
		
		logger.info("GraphLab engine started.");
		return start(mCorePtr);
		
	}

  /**
   * Schedule the execution of an update function on a particular vertex.
   * 
   * @param vertexId
   *          application vertex ID
   * @param updater
   *          updater to execute
   * @throws NullPointerException
   *           if updater was null.
   * @throws NoSuchElementException
   *           if the specified vertex did not exist in the graph that was
   *           passed to {@link #setGraph(Graph)}.
   * @throws IllegalStateException
   *           if {@link #destroy()} was invoked on this object
   */
	public void schedule(int vertexId, Updater updater) {

		if (null == updater)
			throw new NullPointerException("updater must not be null.");
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");

		// map from application vertex ID to graphlab vertex ID
		Integer id = mIdMap.get(vertexId);
		if (null == id)
			throw new NoSuchElementException("vertex did not exist in the graph that was passed to #setGraph.");

		schedule(mCorePtr, updater.ptr(), id);
		
	}
	
	/**
	 * Schedules the execution of the update function on all vertices.
	 * 
	 * @param updater
	 *         updater to execute
   * @throws NullPointerException
   *           if updater was null
   * @throws IllegalStateException
   *           if {@link #destroy()} was invoked on this object
	 */
	public void scheduleAll(Updater updater){
	  
	  if (null == updater)
	    throw new NullPointerException("updater must not be null.");
	  
	  if (mDestroyed)
	    throw new IllegalStateException("Core has been destroyed and may not be reused.");

	  scheduleAll(mCorePtr, updater.ptr());
	  
	}
	
	/**
	 * Sets the number of cpus that the engine will use.
   * This will destroy the current engine and any tasks associated with the current scheduler.
	 * @param ncpus      number of CPUs that the engine will use.
	 */
	public void setNCpus(long ncpus){
	  if (0 >= ncpus)
	    throw new IllegalArgumentException ("ncpus must be positive.");
	  setNCpus(mCorePtr, ncpus);
	}
	
	/**
   * Sets the type of scheduler. 
   * This will destroy the current engine and any tasks currently associated with the scheduler.
   * @param scheduler
   */
  public void setSchedulerType(Scheduler scheduler) {
    setSchedulerType(mCorePtr, scheduler.toString());
  }
	
	/**
   * Sets the scope consistency model that the engine will use. 
   * This will destroy the current engine and any tasks currently associated with the scheduler.
   * @param scope
   */
	public void setScopeType(Scope scope){
	  setScopeType(mCorePtr, scope.toString());
	}
	
	protected Map<Integer, Integer> idMap (){
	  return mIdMap;
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
	 *        {@link #mCorePtr}
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
	 * @param core_ptr
	 *       {@link #mCorePtr}
	 * @param updater_ptr
	 *       {@link org.graphlab.Updater#ptr}
	 * @param vertex_id
	 *       graphlab vertex ID
	 */
	private native void schedule(long core_ptr, long updater_ptr, int vertex_id);

	/**
	 * Add the given function to all vertices using the given priority
	 * @param core_ptr
	 *       {@link #mCorePtr}
   * @param updater_ptr
   *       {@link org.graphlab.Updater#ptr}
	 */
	private native void scheduleAll(long core_ptr, long updater_ptr);
	
  /**
   * Run the engine until a termination condition is reached or there are no
   * more tasks remaining to execute.
   * @param ptr
   *      {@link #mCorePtr}
   * @return runtime
   */
	private native double start(long ptr);

	/**
	 * Set the number of cpus that the engine will use.
	 * This will destroy the current engine and any tasks associated with the current scheduler.
	 * @param ptr
	 *       {@link #mCorePtr}
	 * @param ncpus
	 *       number of CPUs
	 */
	private native void setNCpus(long ptr, long ncpus);
	
	/**
	 * Set the type of scheduler. 
	 * This will destroy the current engine and any tasks currently associated with the scheduler.
	 * @param ptr
	 *       {@link #mCorePtr}
	 * @param schedulerType
	 */
	private native void setSchedulerType(long ptr, String schedulerType);
	
	/**
	 * Set the scope consistency model used in this engine. 
	 * This will destroy the current engine and any tasks associated with the current scheduler.
	 * @param ptr
	 *       {@link #mCorePtr}
	 * @param scopeType
	 */
	private native void setScopeType(long ptr, String scopeType);
	
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
