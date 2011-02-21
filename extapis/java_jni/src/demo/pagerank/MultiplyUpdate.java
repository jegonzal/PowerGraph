/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package demo.pagerank;
 
import graphlab.*;
import graphlab.dense.WeightedEdge;
import graphlab.dense.ScalarVertex;


/**
 * Pagerank-type udpate functions. Computes weighted sum of neighbors'
 * values and and writes its new value to outbound edges.
 * @author akyrola
 */
public class MultiplyUpdate extends UpdateFunction {

    private static final double THRESHOLD = 1E-5;

    public void update(Scope scope, Scheduler scheduler) {
        ScalarVertex vertex = (ScalarVertex) scope.getVertex();
        float newValue = vertex.getSelfEdgeWeight() * vertex.getValue();
        float oldValue = vertex.value;

        /* Calculate new value by multiplying values from incoming
           edges with wege weights
         */
        for(Edge e : scope.getInboundEdges()) {
            WeightedEdge we = (WeightedEdge) e;
            newValue += we.getWeight() * we.getValue();

        }
        vertex.setValue(newValue);

        /**
         * Send my new value to children if it changes significantly
         */
        if (Math.abs(newValue - oldValue) > THRESHOLD) {
             for(Edge outEdge : scope.getOutboundEdges()) {
                 WeightedEdge we = (WeightedEdge) outEdge;
                 we.setValue(newValue);
                 scheduler.addTask(we.getToVertex(), this);
             }
        }
     }
}

