package graphlab.dense;

import graphlab.*;
/**
 * Graph class that has more compact data structure. As a trade-off, it only
 * allows vertices of type DenseVertex.
 *
 * This class is not final quality. TODO.
 * @author akyrola
 */
public class DenseGraph extends SparseGraph {

    public DenseGraph(int n) {
        super(n);
    }

    @Override
    public Vertex addVertex(Vertex vertex) {
        if (!(vertex instanceof DenseVertex))
            throw new IllegalArgumentException("Dense graph vertices must be DenseVertex subclasses!");
        return super.addVertex(vertex);
    }

    @Override
    public Vertex addVertex(int vertexId, Vertex vertex) {
        if (!(vertex instanceof DenseVertex))
            throw new IllegalArgumentException("Dense graph vertices must be DenseVertex subclasses!");
        return super.addVertex(vertexId, vertex);

    }

    public void addWeightedEdge(int from, int to, float weight, float initialValue) {
        int idx = ((DenseVertex) vertices[to]).addIncomingEdge(from, weight, initialValue);
        ((DenseVertex) vertices[from]).addOutgoingEdge(to, idx);
    }

    public void setVertexIncomingData(int to, int indexAtDest, float x) {
        ((DenseVertex) vertices[to]).setIncomingDataAtIndex(indexAtDest, x);
    }

    public void finalizeGraph() {
        for(Vertex v : getVertices()) {
            ((DenseVertex) v).optimizeDataStructures();
        }
    }
}
