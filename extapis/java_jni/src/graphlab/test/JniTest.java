package graphlab.test;

/**
 * Example copied from: http://java.sun.com/developer/onlineTraining/Programming/JDCBook/jniexamp.html#examp
 * @author akyrola
 *         Date: Nov 14, 2010
 */
public class JniTest {

    native int dummy(int i);

    //Load the library
    static {
        System.loadLibrary("graphlabjni");
        System.out.println("Loaded library.");
    }

    public static void main(String args[]) {
    	System.out.println ("According to the JNI library, 1+1=" + new JniTest().dummy(1));
    }
}

