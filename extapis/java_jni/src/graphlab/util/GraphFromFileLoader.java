package graphlab.util;

import graphlab.*;
import graphlab.dense.*;

import java.io.*;

import graphlab.dense.ScalarVertex;

/**
 * @author akyrola
 *         Date: Oct 17, 2009
 */
public class GraphFromFileLoader {

    public static Graph loadStochasticGraph(File f) throws IOException {
        System.out.println("Loading: " + f.getName());
        if (f.getName().endsWith(".csv")) {
            return loadStochasticGraphCSV(f);
        } else {
            return loadStochasticGraphBinary(f);
        }
    }

    private static DenseGraph loadStochasticGraphBinary(File f) throws IOException {
        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 1024000));
        int graphN = dis.readInt();
        long t = System.currentTimeMillis();

        System.out.println("GraphN: " + graphN);

        DenseGraph graph = new DenseGraph(graphN);

        float INITVALUE = 1.0f;

        /* Initialize vertices */
        for(int i=0; i<graphN; i++) {
            graph.addVertex(i, new ScalarVertex(graph, INITVALUE, 500)); // Initialize with uniform distribution
        }

        System.out.println("Vertex init: " + (System.currentTimeMillis() - t) + " ms");
        t = System.currentTimeMillis();

        int i,j;
        float w;
        /* Edges */
        while(true) {
            i = dis.readInt();
            if (i==-1) break;
            j = dis.readInt();
            w = (float) dis.readDouble();
            if (j != i) {
                graph.addWeightedEdge(j, i, w, INITVALUE); // Check, is this right?
            } else {
                ((ScalarVertex) graph.getVertex(i)).setSelfEdgeWeight(w);
            }



        }
        System.out.println("Creating edges: " + (System.currentTimeMillis() - t) + " ms");

        return graph;
    }

    public static DenseGraph loadStochasticGraphCSV(File f) throws IOException  {
        BufferedReader rd = new BufferedReader(new FileReader(f), 128000);

        String ln;
        ln = rd.readLine();
        String[] tok = ln.split(",");

        assert(tok[0].equals("params"));
        int graphN = Integer.parseInt(tok[1]);

        long t = System.currentTimeMillis();

        DenseGraph graph = new DenseGraph(graphN);

        float INITVALUE = 1.0f;

        /* Initialize vertices */
        for(int i=0; i<graphN; i++) {
            graph.addVertex(i, new ScalarVertex(graph, INITVALUE, 100)); // Initialize with uniform distribution
        }

        System.out.println("Vertex init: " + (System.currentTimeMillis() - t) + " ms");
        t = System.currentTimeMillis();

        /* Edges */
        int i, j;
        float weight;
        while((ln = rd.readLine()) != null) {
            tok = ln.split(",");
            i = Integer.parseInt(tok[0]);
            j = Integer.parseInt(tok[1]);
            weight = Float.parseFloat(tok[2]);
            if (j != i) {
                graph.addWeightedEdge(j-1, i-1, weight, INITVALUE); // Check, is this right?
            } else {
                ((ScalarVertex) graph.getVertex(i-1)).setSelfEdgeWeight(weight);
            }
        }
        System.out.println("Creating edges: " + (System.currentTimeMillis() - t) + " ms");

      

        return graph;
    }



}
