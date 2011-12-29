/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package graphlab.dense;

import graphlab.dense.DenseVertex;
import graphlab.dense.DenseGraph;

/**
 * A subclass of DenseVertex carrying a a scalar value. Special handling for self-edges.
 * @author akyrola
 */
public class ScalarVertex extends DenseVertex {
    /** The value associated with the vertex */
    public float value;
    public float selfEdgeWeight;

    /**
     * Create a vertex with some default value
     * @param value
     */
    public ScalarVertex(DenseGraph graph, float value, int estdegree) {
        super(graph, estdegree, estdegree);
        this.selfEdgeWeight = 0.0f;
        this.value = value;
    }

    /**
     * Get the value storted at the vertex
     * @return the value associated with the vertes
     */
    public float getValue() {
        return value;
    }

    /**
     * Set the value of the vertex
     * @param value the new vertex value
     */
    public void setValue(float value){
        this.value = value;
    }

    public float getSelfEdgeWeight() {
        return selfEdgeWeight;
    }

    public void setSelfEdgeWeight(float selfEdgeWeight) {
        this.selfEdgeWeight = selfEdgeWeight;
    }

    public String toString() {
        return value + "";
    }

}
