package graphlab;

/**
 * Simple weighted edge
 * @author akyrola
 */
public class ScalarEdge extends Edge {

    private double weight;

    public ScalarEdge(double weight) {
        this.weight = weight;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }
}
