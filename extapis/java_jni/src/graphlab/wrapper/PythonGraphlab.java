package graphlab.wrapper;

import graphlab.SparseGraph;
import org.python.util.PythonInterpreter;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

/**
 * Application starter for Python/Jython GraphLab applications.
 * Pass python scripts as command-line parameters.
 * <tt>Usage: graphlab.wrapper.PythonGraphlab [loader-script[.py]] [config-script[.pu]][update-function[.py]] [postprocessor[.py]] [inputfile]</tt>
 * @author akyrola
 */
public class PythonGraphlab {

    public static void main(String[] args) {
        if (args.length < 3) {
            System.err.println("Usage: grahlab.wrapper.PythonGraphlab [loader-script[.py]] [config-script[.pu]][update-function[.py]] [postprocessor[.py]] [inputfile]");
        }
        int k = 0;
        String graphloader = args[k++];
        String configScript = args[k++];
        String updateFunction = args[k++];
        String postprocessor = args[k++];

        /* Create empty graph */
        SparseGraph graph = new PythonGraph();

        /* Interpreter */
        long st = System.currentTimeMillis();
        PythonInterpreter interp = new PythonInterpreter();
        interp.set("graph", graph);
        if (args.length > 3) {
            interp.set("filename", args[k++]);
        }

        /* Run graph loader */
        try {
            interp.execfile(new FileInputStream(graphloader));
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + graphloader);
            return;
        }


        System.out.println(":: Graph loaded (" + (System.currentTimeMillis()-st) +  " ms) . ");

        /* Start graphlab */
        PythonUpdateFunctionWrapper wrapfunc = new PythonUpdateFunctionWrapper(updateFunction);

        GraphLabJNIWrapper jniWrapper = new GraphLabJNIWrapper(wrapfunc);

        // TODO: read config
        jniWrapper.setScopeType("vertex");
        jniWrapper.setScheduler("sweep");
        jniWrapper.setMetrics("basic");

        /* Run config */
         try {
            interp.set("graphlab", jniWrapper);
            interp.execfile(new FileInputStream(configScript));
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + configScript);
            return;
        }

        jniWrapper.setGraph(graph);
        jniWrapper.scheduleAll();
        jniWrapper.start();



        /* Run post processor */
        try {
            interp.execfile(new FileInputStream(postprocessor));
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + postprocessor);
            return;
        }
    }

}
