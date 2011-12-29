package demo.scheduling;

import graphlab.*;
import graphlab.wrapper.GraphLabJNIWrapper;

import javax.swing.*;
import java.awt.*;

/**
 * This demo makes a grid and colors the graph randomly.
 * Then each vertex is scheduler to choose a new color, and
 * after update, they schedule their neighbors with random priority.
 * This demo also shows how to use multiple update functions. Other update function
 * changes the color of a pixel to a shade of red, while the other sets it to
 * white.
 * Effect of priority scheduling is evident from the visualization.
 * @author akyrola
 *         Date: Dec 7, 2010
 * TODO update to v2
 */
public class PriorityDemo {

    JFrame frame;
    VisualizerPanel visualizer1;
    SparseGraph graph;

    public PriorityDemo() {
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


    /** Update function 1 **/
    static class ColorUpdateFunction1 extends UpdateFunction {
        private UpdateFunction otherFunction;
        public void update(Scope scope, Scheduler scheduler) {
            ColorVertex colorVertex = ((ColorVertex) scope.getVertex());
            colorVertex.setValue((colorVertex.getValue() + 29)%255);
            // Schedule outbound neighbors
            for(Edge e : scope.getOutboundEdges()) {
                scheduler.addTask(e.getToVertex(), this, (float) Math.random());
            }
            // Schedule inbound neighbors
            for(Edge e : scope.getInboundEdges()) {
                scheduler.addTask(e.getToVertex(), otherFunction, (float) Math.random());
            }
        }

        public void setOtherFunction(UpdateFunction otherFunction) {
            this.otherFunction = otherFunction;
        }
    }

    /** Update function 2 **/
     static class ColorUpdateFunction2 extends UpdateFunction {
        private UpdateFunction otherFunction;
        public void update(Scope scope, Scheduler scheduler) {
            ColorVertex colorVertex = ((ColorVertex) scope.getVertex());
            colorVertex.setValue(0);
            // Schedule outbound neighbors
            for(Edge e : scope.getOutboundEdges()) {
                scheduler.addTask(e.getToVertex(), this, (float) Math.random());
            }
            // Schedule inbound neighbors
            for(Edge e : scope.getInboundEdges()) {
                scheduler.addTask(e.getToVertex(), otherFunction, (float) Math.random());
            }
        }

        public void setOtherFunction(UpdateFunction otherFunction) {
            this.otherFunction = otherFunction;
        }
    }

    /**
     * Setup GraphLab and go.
     */
    void start() {
        ColorUpdateFunction1 updateFunction1 = new ColorUpdateFunction1();
        ColorUpdateFunction2 updateFunction2 = new ColorUpdateFunction2();
        updateFunction1.setOtherFunction(updateFunction2);
        updateFunction2.setOtherFunction(updateFunction1);

        GraphLabJNIWrapper graphlabJni = new GraphLabJNIWrapper(updateFunction1, updateFunction2);
        graphlabJni.setGraph(graph);
        graphlabJni.setScheduler("priority");
        for (int vid=0; vid<graph.getNumOfVertices(); vid++)
            graphlabJni.addTask(vid, updateFunction2);
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
        final PriorityDemo demo = new PriorityDemo();
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
