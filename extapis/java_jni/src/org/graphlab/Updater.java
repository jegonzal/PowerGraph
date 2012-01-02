package org.graphlab;

/**
 * Updater
 * 
 * The GraphLab engine will invoke an updater on each scheduled node. Extend
 * this class to provide an update function for that node. Note that the update
 * function may update node data, modify edge data, and schedule neighbors, but
 * may not modify the graph structure. You may reuse the updater object on
 * across multiple vertices (this is encouraged). Engine uses {@link #id()} to
 * uniquely identify updaters.
 * 
 * @author Jiunn Haur Lim
 */
public abstract class Updater {

	/** Flag to indicate that the updater has no ID */
	public static final int ID_NOT_SET = -1;

	/**
	 * Identifier for this updater Will be set when this updater is first
	 * scheduled.
	 */
	private int mId = ID_NOT_SET;

	/**
	 * Updates the vertex identified by <tt>vertex_id</tt>. Subclasses may wish
	 * to maintain a reference to the graph object.
	 * 
	 * @param contextPtr   address of C++ context object
	 * @param vertexId     application vertex ID
	 */
	public abstract void update(long corePtr, long contextPtr, int vertexId);

	/**
	 * Do not call. This is only useful to the scheduler.
	 * 
	 * @return ID of this updater
	 */
	protected final int id() {
		return mId;
	}

	/**
	 * Do not call. This is only useful to the scheduler.
	 * 
	 * @param id
	 *            ID of this updater
	 */
	protected final void setId(int id) {
		mId = id;
	}

	/* 
	 * Two updaters are equal if they share the same id.
	 * (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals (Object other) {
		if (null == other || !(other instanceof Updater))
			return false;
		return mId == ((Updater) other).id();
	}

}
