package graphlab;

import java.util.List;

/**
 * Scope interface. Instance of Scope is passed to update function, and allows
 * querying vertex neighborhood.
 * @author akyrola
 */
public interface Scope {
    
    /**
      * Get current vertex id.
      */
    public int getVertexId();
    
    /**
      * Get current vertex object.
      */
    Vertex getVertex();
    
    /**
      * @param list of inbound edges (unmodifiable).
      */
    public List<Edge> getInboundEdges();

    /**
      * @param list of out bound edges (unmodifiable).
      */
    public List<Edge> getOutboundEdges();
 
    
    public Vertex getNeighbor(int vertexId);

    public int getNumInboundEdges();

    public int getNumOutboundEdges();

    public int getNumOfVertices();

}
