package graphlab.wrapper;

import graphlab.SparseGraph;
import graphlab.Vertex;
import graphlab.Edge;
import org.python.core.PyObject;

/**
 * Wrapper object for python graphs. Allows storing any python (Jython)
 * type in vertex and edges.
 * @author akyrola
 */
public class PythonGraph extends SparseGraph {


    public int addVertex(PyObject val) {
        return super.addVertex(new PythonVertex(val)).getId();
    }

    public void addEdge(int from, int to, PyObject value) {
        super.addEdge(new PythonEdge(value), from, to);
    }

    

    public static class PythonVertex extends Vertex {
        public PyObject value;

        PythonVertex(PyObject value) {
            this.value = value;
        }
    }

    public static class PythonEdge extends Edge {
        public PyObject value;

        PythonEdge(PyObject value) {
            this.value = value;
        }
    }

}
