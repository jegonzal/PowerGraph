package graphlab.test;

/**
 * Example copied from: http://java.sun.com/developer/onlineTraining/Programming/JDCBook/jniexamp.html#examp
 * @author akyrola
 *         Date: Nov 14, 2010
 */
public class JniTest {




    native int dummy(int i);
    int dummyJava(int j) {return j+1;}

    //Load the library
    static {
        System.loadLibrary("jnitest");
        System.out.println("Loaded library.");
    }

    public static void main(String args[]) {


    }
}

