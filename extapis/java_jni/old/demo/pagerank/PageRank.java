/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package demo.pagerank;

import graphlab.*;
import graphlab.wrapper.GraphLabJNIWrapper;
import graphlab.util.GraphFromFileLoader;

import java.io.File;

import graphlab.dense.ScalarVertex;

/**
 * A simple Pagerank-type application.
 * TODO update to v2
 */
public class PageRank {

    /**
     * @param args the command line arguments
     * @throws java.lang.Exception if something goes wrong.
     */
    public static void main(String[] args) throws Exception {
        // Validate Arguments
        if (args.length < 1) {
            throw new RuntimeException(
                    "Usage: PageRank [graph-data-file-csv]");
        }

        // Load graph from data (written my Matlab script)
        Graph graph = GraphFromFileLoader.loadStochasticGraph(
                new File(args[0]));

        // Provide initial task list to GraphLab
        MultiplyUpdate updateFunction = new MultiplyUpdate();

        GraphLabJNIWrapper jniWrapper = new GraphLabJNIWrapper(updateFunction);
        jniWrapper.setScopeType("vertex");
        jniWrapper.setScheduler("sweep");
        jniWrapper.setGraph(graph);
        jniWrapper.scheduleAll();
        jniWrapper.start();

        // Normalize
        double sum = 0.0;
        for(int vid=0; vid<graph.getNumOfVertices(); vid++) {
           sum += ((ScalarVertex) graph.getVertex(vid)).value;
        }

        // Writes the results in Matlab form
        StringBuffer sb = new StringBuffer(10000);
        sb.append("% Please copypaste to Matlab:\n");
        sb.append("V = [");
        int i = 0;
        for (Vertex v : graph.getVertices()) {
            ScalarVertex sv = (ScalarVertex) v;
            sb.append(sv.getValue()/sum + " ");
            if (i++ % 30 == 0) {
                sb.append("...\n");
            }
        }

        sb.append("]';\n sum(V) \n");
        System.out.println(sb);      
    }
}
