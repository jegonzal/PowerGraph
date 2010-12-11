package graphlab.wrapper;

import graphlab.UpdateFunction;
import graphlab.Scope;
import graphlab.Scheduler;
import org.python.core.PyCode;
import org.python.util.PythonInterpreter;

import java.io.FileReader;

/**
 * Wraps a Python (Jython) update function.
 * @author akyrola
 *         Date: Nov 22, 2010
 */
public class PythonUpdateFunctionWrapper extends UpdateFunction {

    private String filename;
    private PyCode compiledCode;

    public PythonUpdateFunctionWrapper(String filename) {
        this.filename = filename;
        compile();
    }

    private void compile() {
        try {
            PythonInterpreter interp = new PythonInterpreter();
            compiledCode = interp.compile(new FileReader(filename));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    static ThreadLocal<PythonInterpreter> threadLocalInterpreter =
            new ThreadLocal<PythonInterpreter>() {
                @Override
                protected PythonInterpreter initialValue() {
                    return new PythonInterpreter();
                }
            };


    public void update(Scope scope, Scheduler scheduler) {
        // Each thread has to have its own interpreter, otherwise
        // namespace gets screwed.
        PythonInterpreter interp = threadLocalInterpreter.get();
        interp.set("scope", scope);
        interp.set("scheduler", scheduler);
        interp.exec(compiledCode);
    }
}
