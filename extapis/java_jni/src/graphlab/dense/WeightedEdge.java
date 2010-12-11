package graphlab.dense;

import graphlab.Edge;

/**
 * Abstract class for weighted edge with a float value.
 * @author akyrola
 *         Date: Dec 21, 2009
 */
public abstract class WeightedEdge extends Edge {

    public abstract float getWeight();
    public abstract float getValue();
    public abstract void setValue(float v);
    public abstract void setWeight(float w);

}
