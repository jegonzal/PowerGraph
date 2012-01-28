package org.graphlab;

import java.util.Set;

import org.apache.log4j.Logger;
import org.graphlab.data.Vertex;
import org.jgrapht.DirectedGraph;

/**
 * GraphLab Core.
 * 
 * <p>
 * This interfaces with the C++ library via
 * <abbr title="Java Native Interface">JNI</abbr> and
 * mirrors <tt>graphlab::core</tt>.
 * </p>
 * 
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public final class Core {

	/** Logger (Java layer only, the C++ layer has its own logger) */
	private static final Logger logger = Logger.getLogger (Core.class);

	/** Indicates if the core has been destroyed and cannot be reused. */
	private boolean mDestroyed = false;

	/** Address of graphlab::core object */
	private long mCorePtr;

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
   * Creates a new GraphLab core.
   * <b>Call {@link #destroy()} when done to free up resources.</b>
   * 
   * @param config      configuration e.g. scheduler, scope
   * @throws CoreException if there was an error creating the core
   */
  public Core(CoreConfiguration config) throws CoreException {

    mCorePtr = createCore(config.toString());
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
	 *            the graph to operate on.
	 * @param <G>
   *          Graph type must implement {@link org.jgrapht.DirectedGraph}. On the
   *          C++ side, a proxy graph is given to the core engine, and updates are
   *          forwarded to the Java updater (which you will provide through
   *          {@link #schedule(Vertex, Updater)}.)
	 * @throws NullPointerException if graph is null
	 * @throws IllegalArgumentException if graph is empty
	 * @throws IllegalStateException if {@link #destroy()} was invoked on this object
	 */
	public <G extends DirectedGraph<V, E>, V extends Vertex, E> void setGraph(G graph) {

		if (null == graph)
			throw new NullPointerException("graph must not be null.");
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");
		
		// inspect vertices
		Set<? extends Vertex> vertices = graph.vertexSet();
		if (null == vertices || 0 == vertices.size())
			throw new IllegalArgumentException("graph must not be empty.");
		
		long startTime = System.currentTimeMillis();

		// add vertices
		for (Vertex vertex : vertices)
		  vertex.setRawId(addVertex(mCorePtr, vertex));

		// add edges
		for (E edge : graph.edgeSet()) {
		  addEdge(mCorePtr, graph.getEdgeSource(edge).rawId(), graph.getEdgeTarget(edge).rawId());
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

	}

  /**
   * Run the engine until a termination condition is reached or there are no
   * more tasks remaining to execute.
   * 
   * @return runtime
   * @throws IllegalStateException if {@link #destroy()} was invoked on this object
   */
	public double start() {
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");
		
		logger.info("GraphLab engine started.");
		return start(mCorePtr);
		
	}
	
	/**
	 * Get the number of updates executed by the engine.
	 * @return update count
	 */
	public long lastUpdateCount(){
	  
	  if (mDestroyed)
      throw new IllegalStateException("Core has been destroyed and may not be reused.");
	  
	  return lastUpdateCount(mCorePtr);
	  
	}

  /**
   * Schedule the execution of an update function on a particular vertex.
   * 
   * @param vertexId
   *          application vertex ID
   * @param updater
   *          updater to execute
   * @throws NullPointerException
   *           if updater or vertex was null.
   * @throws IllegalStateException
   *           if {@link #destroy()} was invoked on this object
   */
	public void schedule(Vertex vertex, Updater<?> updater) {

		if (null == updater || null == vertex)
			throw new NullPointerException("updater and vertex must not be null.");
		
		if (mDestroyed)
			throw new IllegalStateException("Core has been destroyed and may not be reused.");

		schedule(mCorePtr, updater, vertex.rawId());
		
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
	public void scheduleAll(Updater<?> updater){
	  
	  if (null == updater)
	    throw new NullPointerException("updater must not be null.");
	  
	  if (mDestroyed)
	    throw new IllegalStateException("Core has been destroyed and may not be reused.");

	  scheduleAll(mCorePtr, updater);
	  
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
   * Sets the type of scheduler. This only sets the type, and ignores any
   * scheduler options. This will destroy the current engine and any tasks
   * currently associated with the scheduler.
   * 
   * @param scheduler
   */
  public void setSchedulerType(Scheduler scheduler) {
    setSchedulerType(mCorePtr, scheduler.type());
  }
	
	/**
   * Sets the scope consistency model that the engine will use. 
   * This will destroy the current engine and any tasks currently associated with the scheduler.
   * @param scope
   */
	public void setScopeType(Scope scope){
	  setScopeType(mCorePtr, scope.toString());
	}
	
	/**
	 * Creates and initializes graphlab::core -> dynamically allocates a core.
	 * Must be freed by a corresponding call to {@link #destroyCore()}.
	 * 
	 * @return address of core or 0 on failure
	 */
	private native long createCore();

	/**
   * Creates and initializes graphlab::core -> dynamically allocates a core.
   * Must be freed by a corresponding call to {@link #destroyCore()}. 
   * 
   * @return address of core or 0 on failure
   */
	private native long createCore(String command_line_args);

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
	 * @param vertex
	 * 			application vertex ID
	 * @return
	 * 			graphlab vertex ID
	 */
	private native int addVertex(long ptr, Vertex vertex);

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
	 * @param updater
	 * @param vertex
	 *       graphlab vertex ID
	 */
	private native void schedule(long core_ptr, Updater<?> updater, int vertexId);

	/**
	 * Add the given function to all vertices using the given priority
	 * @param core_ptr
	 *       {@link #mCorePtr}
   * @param updater
	 */
	private native void scheduleAll(long core_ptr, Updater<?> updater);
	
  /**
   * Run the engine until a termination condition is reached or there are no
   * more tasks remaining to execute.
   * @param ptr
   *      {@link #mCorePtr}
   * @return runtime
   */
	private native double start(long ptr);
	
	/**
	 * Gets the number of updates executed by the engine.
	 * @param ptr {@link Core#mCorePtr}
	 * @return update count
	 */
	private native long lastUpdateCount(long ptr);

	/**
	 * Set the number of cpus that the engine will use.
	 * This will destroy the current engine and any tasks associated with the current scheduler.
	 * If this is not what you want, then configure the core in the constructor instead.
	 * 
	 * @param ptr
	 *       {@link #mCorePtr}
	 * @param ncpus
	 *       number of CPUs
	 */
	private native void setNCpus(long ptr, long ncpus);
	
	/**
	 * Set the type of scheduler. 
	 * This will destroy the current engine and any tasks currently associated with the scheduler.
	 * If this is not what you want, then configure the core in the constructor instead.
	 * 
	 * @param ptr
	 *       {@link #mCorePtr}
	 * @param schedulerType
	 */
	private native void setSchedulerType(long ptr, String schedulerType);
	
	/**
	 * Set the scope consistency model used in this engine. 
	 * This will destroy the current engine and any tasks associated with the current scheduler.
	 * If this is not what you want, then configure the core in the constructor instead.
	 * 
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
