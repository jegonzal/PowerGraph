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
      * @return list of inbound edges (unmodifiable).
      */
    public List<Edge> getInboundEdges();

    /**
      * @return list of out bound edges (unmodifiable).
      */
    public List<Edge> getOutboundEdges();
 
    
    /**
      * @return neighbor vertex object
      */
    public Vertex getNeighbor(int vertexId);

    public int getNumInboundEdges();

    public int getNumOutboundEdges();
    
    /**
      * @return number of vertices in graph
      */
    public int getNumOfVertices();

}
