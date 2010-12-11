package graphlab;


/**
 * Efficient representation for sparse graphs. Default
 * Graph class to be used.
 * @author akyrola
 */
public class SparseGraph implements Graph {

    protected Vertex[] vertices;

    int n = 0, maxId = 0;
    boolean coloringenabled = false;

    /**
     *  Constructor.
     */
    public SparseGraph() {
        this(1000);
    }


    /**
     * Creates a new graph with estimated number of vertices.
     *  Internal data structures will grow if the estimate is exceeded.
     * @param n estimated number of vertices
     */
    public SparseGraph(int n) {
        vertices = new Vertex[n];
    }


    /**
     * Returns Vertex object of given id.
     * @param vid vertex id (>=0)
     * @return vertex object
     * @throws java.lang.ArrayIndexOutOfBoundsException
     */
    public Vertex getVertex(int vid) {
        if (vid<0 || vid>=vertices.length) throw new ArrayIndexOutOfBoundsException("Vertex id out of bounds!");
        return vertices[vid];
    }


    /**
     * Add a vertex to the graph.
     * @param vertex
     * @return the vertex object. Now id-field is set
     */
    public Vertex addVertex(Vertex vertex) {
        return addVertex(n, vertex);
    }

    protected synchronized Vertex addVertex(int vertexId, Vertex vertex) {
        checksize(vertexId);
        if (vertices[vertexId] != null) throw new IllegalArgumentException("Tried to change vertex data! Vertex=" + vertexId);
        vertices[vertexId] = vertex;
        vertex.setId(vertexId);
        n++;
        maxId = Math.max(vertexId, maxId);
        return vertex;
    }

    /**
     * Set vertex color. Graph coloring can be used for scheduling with the "colored_scheduler"
     * scheduler.
     * @param vertex vertex id
     * @param vertex_color color id. Color must be integer between 0 and 255.
     */
    public void setVertexColor(int vertex, int vertex_color) {
        if (vertex_color >= 256) throw new IllegalArgumentException("Maximum vertex color is 255!");
        getVertex(vertex).vertexcolor = (short) vertex_color;
        coloringenabled = true;
    }

    /**
     * Add a directed  edge
     * @param e edge object encapsulating edge data
     * @param fromVertex vertex id of the source vertex
     * @param toVertex vertex id of the target vertex
     */
    public synchronized void addEdge(Edge e, int fromVertex, int toVertex) {
        if (fromVertex >= vertices.length || vertices[fromVertex] == null)
            throw new IllegalArgumentException("Tried to add edge to non-existing vertex " + fromVertex);
        if (toVertex >= vertices.length || vertices[toVertex] == null)
            throw new IllegalArgumentException("Tried to add edge to non-existing vertex " + toVertex);
        e.setFromVertex(fromVertex);
        e.setToVertex(toVertex);

        vertices[fromVertex].addOutboundEdge(e);
        vertices[toVertex].addInboundEdge(e);
    }


    /**
     * Returns a <b>copy</b> of vertex ids of a graph.
     * @return
     */
    public int[] vertices() {
        int[] v = new int[n];
        int j = 0;
        for(int i=0; i< vertices.length; i++) {
            if (vertices[i] != null) v[j++] = i;
        }
        return v;
    }

    /**
     * Returns list of outbound neighbors of a vertex. Do not modify.
     * @param vertex  id of vertex
     */
    public int[] children(int vertex) {
        return vertices[vertex].getChildren();
    }


    /**
     * Returns list of inbound neighbors of a vertex. Do not modify.
     * @param vertex  id of vertex
     */
    public int[] parents(int vertex) {
        return vertices[vertex].getParents();
    }


    /**
     * Grows the data structures if needed
     * @param vertexId
     */
    protected synchronized void checksize(int vertexId) {
        if (n<vertexId) n = vertexId;
        if (vertexId >= vertices.length) {
            /* Grow by 10% */
            int newCapacity = (int) (vertexId * 1.1);
            Vertex[] newVertexData = new Vertex[newCapacity];
            System.arraycopy(vertices, 0, newVertexData, 0, vertices.length);

            vertices = newVertexData;
        }
    }


    /**
     * Returns a <b>copy</b> of vertex-array of the graph.
     * @return array of vertex objects.
     */
    public Vertex[] getVertices() {
          Vertex[] vcopy = new Vertex[maxId+1];
          System.arraycopy(vertices, 0, vcopy, 0, maxId+1);
          return vcopy;
    }


    // TODO: sort edges and use binary search
    /**
     * Returns an edge from source vertex to target vertex. Note:
     * this is a slow operation (O(degree)) since edges are not sorted.
     * You are usually best of using iterators from getVertex(vid).getInboundEdges()
     *
     * @param sourceVertexId  vertex id of the origin vertex
     * @param targetVertexId  vertex id of the destination vertex
     * @return edge spanning source to target vertex
     */
    public Edge getEdge(int sourceVertexId, int targetVertexId) {
        return vertices[sourceVertexId].getEdgeTo(targetVertexId);
    }


    /**
      * Whether this graph is a "colored" graph. Graph coloring
      * can be used for scheduling.
      */
    public boolean isColoringEnabled() {
        return coloringenabled;
    }


    /**
     * Returns the number of vertices in graph.
     * @return
     */
    public int getNumOfVertices() {
        return n;
    }

    public void finalizeGraph() {
        // not needed.
    }

    public String toString() {
        StringBuffer sb = new StringBuffer(vertices.length * 30);
        for(Vertex v : vertices) {
            if (v != null) {
                sb.append(v.getId());
                sb.append(", ");
                sb.append(v);
                sb.append(", ");
                sb.append(v.getOutboundEdges());
                sb.append("\n");
            }
        }
        return sb.toString();
    }
}
