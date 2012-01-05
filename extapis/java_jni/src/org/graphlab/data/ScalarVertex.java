package org.graphlab.data;

/**
 * Vertex carrying a single scalar value. Special handling for self-edges by
 * storing self-edge weights.
 * 
 * @author jegonzal
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class ScalarVertex implements Vertex {

	/** Value associated with the vertex */
	private float mValue;

	/** Weight of a possible self-edge */
	private float mSelfEdgeWeight = 0.0f;

	/** Vertex ID */
	private int mId = 0;

	/**
	 * Create a vertex with some default value.
	 * 
	 * @param value
	 *            default value
	 */
	public ScalarVertex(float value) {
		mValue = value;
	}

	/**
	 * Gets the value stored at the vertex
	 * 
	 * @return the value associated with the vertex
	 */
	public float value() {
		return mValue;
	}

	/**
	 * Sets the value of the vertex
	 * 
	 * @param value
	 *            the new vertex value
	 */
	public void setValue(float value) {
		mValue = value;
	}

	/**
	 * Gets the weight of a possible self-edge
	 * 
	 * @return weight
	 */
	public float selfEdgeWeight() {
		return mSelfEdgeWeight;
	}

	/**
	 * Sets the weight of a possible self-edge
	 * 
	 * @param weight
	 */
	public void setSelfEdgeWeight(float weight) {
		mSelfEdgeWeight = weight;
	}

	@Override
	public String toString() {
		return mValue + "";
	}

	/**
	 * @return vertex id
	 */
	public int id() {
		return mId;
	}

	/**
	 * Sets the vertex id
	 * 
	 * @param id
	 *            vertex id
	 */
	public void setId(int id) {
		mId = id;
	}

	@Override
	public boolean equals(Object other) {
		if (null == other || !(other instanceof ScalarVertex))
			return false;
		return mId == ((ScalarVertex) other).id();
	}

}
