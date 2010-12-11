/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package graphlab;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Vertex class. Subclass to encapsulate your application specific data.
 * @author akyrola
 */
public class Vertex {

    private LinkedList<Edge> inboundEdges  = new LinkedList<Edge>();
    private LinkedList<Edge> outboundEdges = new LinkedList<Edge>();
    protected int id;
    short vertexcolor = 0;

    /**
     * @return unmodifiable list of inbound edges.
     */
    public List<Edge> getInboundEdges() {
        return Collections.unmodifiableList(inboundEdges);
    }

    /**
     *
     * @return unmodifiable list of outbound edges
     */
    public List<Edge> getOutboundEdges() {
        return Collections.unmodifiableList(outboundEdges);
    }

    /**
     * Use Graph.addEdge()
     * @param e
     */
    void addInboundEdge(Edge e) {
        inboundEdges.add(e);
    }


    /**
     * Use Graph.addEdge()
     * @param e
     */
    void addOutboundEdge(Edge e) {
        outboundEdges.add(e);
    }

    /**
     * Returns edge directed to vertex toVertexId. Note: slow operation.
     * @param toVertexId
     */
    public Edge getEdgeTo(int toVertexId) {
        for(Edge e:outboundEdges) {
            if (e.to == toVertexId) {
                return e;
            }
        }
        return null;
    }

    public Edge getEdgeFrom(int fromVertexId) {
        for(Edge e:inboundEdges) {
            if (e.from == fromVertexId) {
                return e;
            }
        }
        return null;
    }

    /**
     * @return list of child vertices (targets of outbound edges)
     */
    public int[] getChildren() {
        // TODO: optimize by caching the list
        int[] v = new int[outboundEdges.size()];
        int i = 0;
        for(Edge e : outboundEdges) v[i++] = e.getToVertex();
        return v;
    }


    /**
     * @return list of parent vertices (targets of inbound edges)
     */
    public int[] getParents() {
        // TODO: optimize by caching the list
        int[] v = new int[inboundEdges.size()];
        int i = 0;
        for(Edge e : inboundEdges) v[i++] = e.getFromVertex();
        return v;
    }

    /**
     * Returns vertex id
     * @return
     */
    public int getId() {
        return id;
    }


    void setId(int id) {
        this.id = id;
    }


    /**
     * Returns vertex color.
     * @return
     */
    public int getVertexColor() {
        return vertexcolor;
    }





} // End Vertex
