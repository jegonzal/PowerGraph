package demo.scheduling;

import graphlab.*;
import graphlab.wrapper.GraphLabJNIWrapper;

import javax.swing.*;
import java.awt.*;

/**
 * This demo makes a grid and colors the graph so that adjacent vertices
 * have different color (checker board coloring). Then a colored scheduler
 * is run, which updates vertices with same color in parallel. Each vertex
 * picks a value randomly, and a visualization shows how the computation proceeds.
 * @author akyrola
 *         Date: Dec 4, 2010
 * TODO update to v2
 */
public class ColoredSchedulerDemo {

    JFrame frame;
    VisualizerPanel visualizer1;
    SparseGraph graph;

    public ColoredSchedulerDemo() {
        openWindow();
    }


    /* =========== GRAPHLAB STUFF =============== */
    void initGraph(int width, int height) {
        int numPixels = width * height;
        graph = new SparseGraph(numPixels);

        // Create grid vertices
        for(int i=0; i<numPixels; i++) {
            int row = i/width;
            int col = i%width;
            int color = (row+col)%2;
            graph.addVertex(new ColorVertex(120*color));
            graph.setVertexColor(i, color);
        }
        // Connect vertices with color 0 to vertices with color
         for(int i=0; i<numPixels; i++) {
            int row = i/width;
            int col = i%width;
            ColorVertex colorVertex = ((ColorVertex) graph.getVertex(i));
            if (colorVertex.getVertexColor() ==0) {
                if (row > 0) graph.addEdge(new ScalarEdge(1.0), (row-1)*width + col, i);
                if (col > 0) graph.addEdge(new ScalarEdge(1.0), i-1, i);
                if (col+1 < width) graph.addEdge(new ScalarEdge(1.0), i+1, i);
                if (row+1 < height) graph.addEdge(new ScalarEdge(1.0), (row+1)*width + col, i);
            }
        }
    }

    static class ColorUpdateFunction extends UpdateFunction {
        public void update(Scope scope, Scheduler scheduler) {
           ColorVertex colorVertex = ((ColorVertex) scope.getVertex());
            // Set a new value for vertices of different color with a bit
            // different way. This way it is easier to see what is happening.
           switch(colorVertex.getVertexColor()) {
               case 0:
                   // Take average
                   int sum = 0;
                   for(Edge e : scope.getInboundEdges()) {
                       sum += ((ColorVertex) scope.getNeighbor(e.getFromVertex())).getValue();
                   }
                   colorVertex.setValue(sum / scope.getInboundEdges().size());
                   break;
               case 1:
                   colorVertex.setValue((colorVertex.getValue() + 29)%255);
           }
        }
    }

    /**
     * Setup GraphLab and go.
     */
    void start() {
        GraphLabJNIWrapper graphlabJni = new GraphLabJNIWrapper(new ColorUpdateFunction());
        graphlabJni.setGraph(graph);
        graphlabJni.setScheduler("colored");
        graphlabJni.scheduleAll();
        graphlabJni.setIterations(5000);
        graphlabJni.start();
    }


    /* =========== GUI / WINDOW STUFF =========== */
    private void openWindow() {
        frame = new JFrame("GraphLab visualizer client");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        frame.pack();
        frame.setVisible(true);
    }



    void initializeViews() {
        int num_vertices = graph.getNumOfVertices();
        int width = (int) Math.sqrt(num_vertices);
        int height = num_vertices/width+1;


        /* Auto zoom */
        float wzoom = (320.0f/width);
        float hzoom = (500.0f/height);
        float zoom = Math.min(wzoom, hzoom);

        /* Create two visualizer views */
        frame.setSize((int) (width * zoom *3), (int) (height * zoom));

        System.out.println("Num vertices: " + num_vertices + " width: " + width + " height: " + height);
        visualizer1 = new VisualizerPanel(width, height, (int) (zoom*width), (int) (zoom*height) + 50);

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.LINE_AXIS));

        frame.getContentPane().add(topPanel, BorderLayout.CENTER);
        topPanel.add(visualizer1, BorderLayout.CENTER);
        frame.pack();
    }


    void update() {
        visualizer1.updateData(graph);
    }

    public static void main(String[] args) {
        final ColoredSchedulerDemo demo = new ColoredSchedulerDemo();
        demo.initGraph(200,200);
        demo.initializeViews();
        demo.update();

        // Start a thread that updates the visuals
        Thread updateThread = new Thread() {
            public void run() {
                long t = System.currentTimeMillis();
                while(System.currentTimeMillis()-t < 50000) {
                    demo.update();
                    try {
                        Thread.sleep(5);
                    } catch (InterruptedException e) {

                    }
                }
            }

        };
        updateThread.start();

        demo.start();
        
    }

}
