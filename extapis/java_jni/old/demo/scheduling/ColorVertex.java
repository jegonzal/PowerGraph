package demo.scheduling;

import graphlab.Vertex;

/**
 * Simple vertex carrying and integer value.
 * @author akyrola
 *         Date: Dec 4, 2010
 * TODO update to v2
 */
public class ColorVertex extends Vertex {

    int value;

    public ColorVertex(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public void setValue(int value) {
        this.value = value;
    }
}
