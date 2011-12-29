package org.graphlab;

public abstract class Updater {

	public static final int ID_NOT_SET = -1;
	
	/** Identifier */
	private int mId = ID_NOT_SET;
	
	/**
	 * Updates the vertex identified by <tt>vertex_id</tt>.
	 * @param vertex_id
	 */
	public abstract void update (int vertex_id);
	
	/**
	 * Do not call.
	 * @return id of this updater
	 */
	protected int id(){
		return mId;
	}
	
	/**
	 * Do not call.
	 * @param id 		ID of this updater
	 */
	protected void setId(int id){
		mId = id;
	}
	
	@Override
	public boolean equals (Object other){
		if (null == other || !(other instanceof Updater)) return false;
		return mId == ((Updater) other).id();
	}
	
}
