package graphlab.test;

/**
 * Example copied from: http://java.sun.com/developer/onlineTraining/Programming/JDCBook/jniexamp.html#examp
 * @author akyrola
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
public class JniTest {

    native int dummy(int i);

    // load the library
    static {
    	System.out.println("Library path: " + System.getProperty("java.library.path"));
        System.loadLibrary("graphlabjni");
        System.out.println("Loaded library.");
    }

    public static void main(String args[]) {
    	System.out.println ("According to the JNI library, 1+1=" + new JniTest().dummy(1));
    }
    
}

