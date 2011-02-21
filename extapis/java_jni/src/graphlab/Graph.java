/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package graphlab;


/**
 * Graph interface.
 * GraphLab graphs are directional. Self-edges are not allowed.
 */
public interface Graph {

  
    /**
     * @return array of vertices
     */
    Vertex[] getVertices();

    /**
     * Returns a vertex id
     * @param vid
     * @return
     */
    Vertex getVertex(int vid);

    /**
     * Add directed edge between two vertices.
     * @param e edge object encapsulating the edge data
     * @param fromVertex
     * @param toVertex
     */
    void addEdge(Edge e, int fromVertex, int toVertex);

     /**
     * Get an edge object.
     * @param sourceVertexId  vertex id of the origin vertex
     * @param targetVertexId  vertex id of the destination vertex
     * @return edge-object corresponding to the edge or null.
     */
    Edge getEdge(int sourceVertexId, int targetVertexId);

    /**
     * Add a vertex to the graph.
     * @param vertex object containing vertex data
     * @return vertex object of the new vertex
     */
    Vertex addVertex(Vertex vertex);


    /**
     * Set vertex color.
     * @param vertex
     * @param vertex_color
     */
    void setVertexColor(int vertex, int vertex_color);


    /**
     * @param vertex  id of vertex
     * @return  array of outbound neighbors of a vertex
     */
    int[] children(int vertex);

    /**
     * @param vertex  id of vertex
     * @return  array of inbound neighbors of a vertex
     */
    int[] parents(int vertex);


    /**
     * Query whether graph is colored.
     * @return true if vertices have color
     */
    boolean isColoringEnabled();


    /**
     * Return the number of vertices in graph.
     * @return
     */
    int getNumOfVertices();

    void finalizeGraph();


} // End of Graph
