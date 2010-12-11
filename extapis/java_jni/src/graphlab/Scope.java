package graphlab;

import java.util.List;

/**
 * @author akyrola
 */
public interface Scope {

    public int getVertexId();

    public List<Edge> getInboundEdges();

    public List<Edge> getOutboundEdges();

    public int[] getOutboundEdgeIds();

    public int[] getInboundEdgeIds();


    public Vertex getNeighbor(int vertexId);

    public int getNumInboundEdges();

    public int getNumOutboundEdges();

    public int getNumOfVertices();

    Vertex getVertex();
}
