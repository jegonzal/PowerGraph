package graphlab.dense;

import graphlab.Vertex;
import graphlab.Edge;

import java.util.List;
import java.util.ArrayList;

/**
 * Vertex with more compact data structure. Allows only float data in edges.
 * Edge-objects are created on-the-fly.
 * @author akyrola
 *         Date: Dec 21, 2009
 */
public class DenseVertex extends Vertex {

    private int inDegree, outDegree;

    private float[] inedgeData;
    private float[] inedgeWeight;
    private int[] inedgeVertexId;
    private int[] outedgeVertexIds;
    private int[] outedgeValueIndices; // index of inedgeData[] in destination

    private final DenseGraph graph;

    public DenseVertex(DenseGraph graph, int estimatedInDegree, int estimatedOutDegree) {
        inedgeData = new float[estimatedInDegree];
        inedgeWeight = new float[estimatedInDegree];
        inedgeVertexId = new int[estimatedInDegree];
        outedgeVertexIds = new int[estimatedOutDegree];
        outedgeValueIndices = new int[estimatedOutDegree];
        this.graph = graph;
    }


    /**
     * Adds an incoming edge with given weight and initial data
     * @param weight
     * @param initialData
     * @return index of new edge's data in the data array of this vertex
     */
    public int addIncomingEdge(int vertexId, float weight, float initialData) {
        int index = inDegree;
        inDegree++;
        checkIndegreeStorage();

        inedgeData[index] = initialData;
        inedgeWeight[index] = weight;
        inedgeVertexId[index] = vertexId;

        return index;
    }

    void setIncomingDataAtIndex(int index, float newData) {
        inedgeData[index] = newData;
    }

    public void addOutgoingEdge(int vertexId, int indexAtDest) {
        int index = outDegree;
        outDegree++;
        checkOutDegreeStorage();

        outedgeVertexIds[index] = vertexId;
        outedgeValueIndices[index] = indexAtDest;

    }


    private void checkIndegreeStorage() {
        if (inDegree >= inedgeData.length) {
            int newSize = (int) (inDegree * 1.2);
            resizeInArrays(newSize);
        }
    }

    private void resizeInArrays(int newSize) {
        if (newSize == inedgeData.length) return;
        float[] tmpdata = inedgeData;
        float[] tmpweights = inedgeWeight;
        int[] tmpvertexids = inedgeVertexId;
        inedgeData = new float[newSize];
        int len = (newSize > tmpdata.length ? tmpdata.length : newSize);
        System.arraycopy(tmpdata, 0, inedgeData, 0, len);
        inedgeWeight = new float[newSize];
        System.arraycopy(tmpweights, 0, inedgeWeight, 0, len);
        inedgeVertexId = new int[newSize];
        System.arraycopy(tmpvertexids, 0, inedgeVertexId, 0, len);
    }

    private void checkOutDegreeStorage() {
        if (outDegree >= outedgeValueIndices.length) {
            int newSize = (int) (outDegree * 1.2);
            resizeOutArrays(newSize);
        }
    }

    private void resizeOutArrays(int newSize) {
        if (newSize == outedgeValueIndices.length) return;
        int[] tmpVertexIds = outedgeVertexIds;
        int[] tmpValueIndices = outedgeValueIndices;
        int len = (newSize > tmpVertexIds.length ? tmpVertexIds.length : newSize);
        outedgeVertexIds = new int[newSize];
        outedgeValueIndices = new int[newSize];
        System.arraycopy(tmpVertexIds, 0, outedgeVertexIds, 0, len);
        System.arraycopy(tmpValueIndices, 0, outedgeValueIndices, 0, len);
    }


    /* Shrinks all internal arrays to proper size */
    public void optimizeDataStructures() {
        resizeInArrays(inDegree);
        resizeOutArrays(outDegree);

    }


    class TempOutEdge extends WeightedEdge {

        int indexAtDest;
        TempOutEdge(int from, int destVertexId, int indexAtDest) {
            this.from = from;
            this.to = destVertexId;
            this.indexAtDest = indexAtDest;
         }

        public void setValue(float x) {
            graph.setVertexIncomingData(to, indexAtDest, x);
        }

        @Override
        public float getWeight() {
            throw new UnsupportedOperationException();
        }

        @Override
        public float getValue() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void setWeight(float w) {
            throw new UnsupportedOperationException();
        }
    }

    class TempInEdge extends WeightedEdge {

        private int idx;

        TempInEdge(int from, int to, int idx) {
            this.from = from;
            this.to = to;
            this.idx = idx;
        }

        public float getWeight() {
            return inedgeWeight[idx];
        }

        public float getValue() {
            return inedgeData[idx];
        }

        @Override
        public void setValue(float v) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void setWeight(float w) {
            throw new UnsupportedOperationException();
        }
    }


    @Override
    public List<Edge> getInboundEdges() {
        ArrayList<Edge> l = new ArrayList<Edge>(inDegree);
        // TODO: thread-local object pooling!
        for(int i=0; i<inDegree; i++) l.add(new TempInEdge(inedgeVertexId[i], this.id, i));
        return l;
    }

    @Override
    public List<Edge> getOutboundEdges() {
        // TODO: thread-local object pooling!        
        ArrayList<Edge> l = new ArrayList<Edge>(outDegree);
        for(int i=0; i<outDegree; i++) l.add(new TempOutEdge(this.id, outedgeVertexIds[i], outedgeValueIndices[i]));
        return l;
    }


    @Override
    public Edge getEdgeTo(int toVertexId) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Edge getEdgeFrom(int fromVertexId) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int[] getChildren() {
        return outedgeVertexIds;
    }

    @Override
    public int[] getParents() {
        throw new UnsupportedOperationException();
    }
}
