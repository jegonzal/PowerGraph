package graphlab.wrapper;

import graphlab.*;

import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Graphlab JNI Wrapper. This will load GraphLab C++ library for multicore and interface with
 * it. Ensure that the JNI library is in the Java's library path (use command line argument -Djava.library.path=...).
 *
 * To modify the number of threads to launch, pass System option as follows: -Dgraphlab.ncpus=8.
 * By default, GraphLab will use all the cores in the system.
 * @author akyrola
 */
public class GraphLabJNIWrapper {



    private AtomicInteger counter = new AtomicInteger(0);

    private List<UpdateFunction> updateFunctions;
    private Graph graph;

    // TODO: support many update functions
    /**
     * Create a Graphlab JNI (Java Native Interface) wrapper.
     * @param updateFunctions  list of update functions that the algorithm will use.
     */
    public GraphLabJNIWrapper(UpdateFunction ...updateFunctions){
        initGraphlab();
        this.updateFunctions = new ArrayList<UpdateFunction>();

        for(int i=0; i<updateFunctions.length; i++) {
            this.updateFunctions.add(updateFunctions[i]);
            updateFunctions[i].functionId = i;
        }

        int ncpus = Integer.parseInt(System.getProperty("graphlab.ncpus",
                Runtime.getRuntime().availableProcessors() + ""));
        setNumCPUs(ncpus);
    }




    /**
     * Set the maximum number of updates to perform.
     * @param taskcount
     */
    native public void setTaskBudget(int taskcount);

    /**
     * Set maximum number of iterations (only valid for some schedulers, such as roundrobin scheduler).
     * @param maxiter
     */
    native public void setIterations(int maxiter);

    /**
     * Set the number of parallel threads to use. Default = all.
     * @param ncpus
     */
    native public void setNumCPUs(int ncpus);

    /**
     * Set the name of the scheduler to use. Examples: multiqueue_fifo, colored, sweep, priority,
     * round_robin, sampling, multiqueue_priority
     * @param s
     */
    native public void setScheduler(String s);
    native public void setScopeType(String s);
    native public void setMetrics(String s);

    /**
     * If using a colored scheduler, use this to automatically compute a
     * greedy graph coloring for the graph. You can also do coloring
     * by hand with graph.setColor().
     */
    native public void computeGraphColoring();



    /**
     * Pass the graph. You must have created the whole graph at this point, and cannot add vertices
     * or edges after this.
     * @param graph
     */
    public void setGraph(Graph graph) {
        long st = System.currentTimeMillis();
        // Initialize graph
        graph.finalizeGraph();
        createGraph(graph.getNumOfVertices());

        // Create edges
        for(int i=0; i<graph.getNumOfVertices(); i++) {
            Vertex v = graph.getVertex(i);
            addEdges(v.getId(), v.getChildren());
        }
        this.graph = graph;
        long t = System.currentTimeMillis()-st;
        System.out.println("Graph transferring took: " + t + " ms.");

        // Set graph coloring
        if (graph.isColoringEnabled()) {
            int[] colors = new int[graph.getNumOfVertices()];
            for(int v=0; v<graph.getNumOfVertices(); v++) {
                colors[v] = graph.getVertex(v).getVertexColor();
            }
            setVertexColors(colors);
        }
    }



    /**
     * Convenience function for adding the first update function to the
     * schedule for all vertices.
     */
    public void scheduleAll() {
        scheduleAll(0);
    }


    /**
     * Start computing. Returns when ready or an error occurs.
     */
    public void start() {
        runGraphlab();
        System.out.println("Update count: " + counter);
    }




    private Scope getScope(final int vertexId) {

        return new Scope() {
            public int getVertexId() {
                return vertexId;
            }

            public List<Edge> getInboundEdges() {
                return graph.getVertex(vertexId).getInboundEdges();
            }

            public List<Edge> getOutboundEdges() {
                return   graph.getVertex(vertexId).getOutboundEdges();
            }

            public int[] getInboundEdgeIds() {
                return graph.parents(vertexId);
            }

            public int[] getOutboundEdgeIds() {
                return graph.children(vertexId);
            }

            // TODO: should I check for neighborness?
            public Vertex getNeighbor(int nbid) {
                return graph.getVertex(nbid);
            }

            public Vertex getVertex() {
                return graph.getVertex(vertexId);
            }

            public int getNumOutboundEdges() {
                return graph.children(vertexId).length;
            }

            public int getNumOfVertices() {
                return graph.getNumOfVertices();
            }

            public int getNumInboundEdges() {
                return graph.parents(vertexId).length;
            }
        };
    }

    class JniScheduler implements Scheduler {

        private int[] vs;
        private int[] funcs;
        private int[] prios;

        int counter = 0;
        UpdateFunction callee;

        JniScheduler(int estsize, UpdateFunction callee) {
            if (estsize<=0) estsize = 1;
            vs = new int[estsize];
            funcs = new int[estsize];
            prios = new int[estsize];
            this.callee = callee;
        }


        public void addTask(int vertex, UpdateFunction function) {
            addTask(vertex, function, 0.0f);
        }

        public void addTask(int vertex, UpdateFunction function, float priority) {
            vs[counter] = vertex;
            funcs[counter] = function.functionId;
            prios[counter] = (int) (1e6 * priority);
            counter++;
            if (counter>=vs.length) {
                int[] tmpvs = new int[counter*2];
                int[] tmpfs = new int[counter*2];
                int[] tmpprios = new int[counter*2];
                System.arraycopy(vs, 0, tmpvs, 0, vs.length);
                System.arraycopy(funcs, 0, tmpfs, 0, funcs.length);
                System.arraycopy(prios, 0, tmpprios, 0, prios.length);
                vs = tmpvs;
                funcs = tmpfs;
                prios = tmpprios;
            }
        }

        // TODO: use System.arraycopy        
        public void addTask(UpdateFunction function, int... v) {
            for(int i : v) addTask(i, function);
        }

        // TODO: use System.arraycopy
        public void addTask(int... v)  {
            addTask(callee, v);
        }

        // TODO: use System.arraycopy
        public void addTask(List<Integer> vv) {
            for(int i : vv) addTask(i, callee);
        }

        // TODO: use System.arraycopy
        public void addTaskToOutbound(Scope scope) {
            addTask(scope.getOutboundEdgeIds());
        }

        int[] toIntArray() {
            int[] arr = new int[counter * 3];
            System.arraycopy(vs, 0, arr, 0, counter);
            System.arraycopy(funcs, 0, arr, counter, counter);
            System.arraycopy(prios, 0, arr, counter*2, counter);

            return arr;
        }
    }


    public void addTask(int vertexid, UpdateFunction func) {
        schedule(new int[]{vertexid}, new int[]{func.functionId});
    }

    /**
     * Internal method. Called by the  GraphLab wrapper in order to execute an update function
     * on a vertex.
     * @param vertexId
     * @param functionId
     * @return new scheduled tasks
     */
    int[] execUpdate(int vertexId, int functionId) {
        try {
            counter.incrementAndGet();
            Scope scope = this.getScope(vertexId);
            JniScheduler scheduler = new JniScheduler(Math.max(1, scope.getNumOutboundEdges()),
                    this.updateFunctions.get(functionId));
            this.updateFunctions.get(functionId).update(scope, scheduler);
            return scheduler.toIntArray();
        } catch (Exception err) {
            err.printStackTrace();
            return new int[0];
        }
    }

    static {
        System.loadLibrary("graphlabjni");
    }

    // Native methods
    native void initGraphlab();
    native void createGraph(int numOfVertices);
    native void addEdges(int fromVertex, int[] toVertex);
    native void runGraphlab();
    native void schedule(int[] vertices, int[] funcIds);
    native void scheduleAll(int funcId);
    native void setVertexColors(int[] colors);



}

