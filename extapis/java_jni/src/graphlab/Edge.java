/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graphlab;


/**
 * This class represents basic directed edge data
 * @author akyrola
 */
public class Edge {

    /** The vertex that the edge originates from */
    public int from;

    /** The vertex that the edge arrives at */
    public int to;

    /** Default constructor  */
    public Edge() {
        from = -1;
        to = -1;
    }


    /**
     * Get the from vertex
     * @return Vertex at the from
     */
    public int getFromVertex() {
        return from;
    }

    /**
     * Get the destination vertex
     * @return Verte at the denstionation
     */
    public int getToVertex() {
        return to;
    }


    void setFromVertex(int from) {
        this.from = from;
    }

    void setToVertex(int to) {
        this.to = to;
    }


    public String toString() {
        return from + " -> " + to;
    }


} // end of Edge



