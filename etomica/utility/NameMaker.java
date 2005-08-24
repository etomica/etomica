/*
 * History
 * Created on Nov 15, 2004 by kofke
 */
package etomica.utility;

import java.util.HashMap;

import etomica.Phase;
import etomica.integrator.Integrator;

/**
 * Utility class with a method that generates a unique string that can be used
 * to name an instance of a class. The name is formed from the class name given
 * to the method, concatenated by a number equal to the number of times the
 * method has been invoked previously with the given class.
 */
public class NameMaker implements java.io.Serializable {

    /**
     * Private constructor prevents instantiation.
     */
    private NameMaker() {
    }

    /**
     * Generates the name String as from the given Class instance. The name is
     * formed from the class name given to the method, concatenated by a number
     * equal to the number of times the method has been invoked previously for
     * the given class.
     */
    public static String makeName(Class className) {
        ClassCount thisCount = (ClassCount) classList.get(className);
        if (thisCount == null) {
            thisCount = new ClassCount();
            classList.put(className, thisCount);
        }
        return className.getName() + thisCount.count++;
    }

    private static HashMap classList = new HashMap();

    //A mutable Integer class used to keep an instance count for each Class.
    private static class ClassCount implements java.io.Serializable {
        int count = 0;
    }

    /**
     * Simple main method to test the function of the class. Should print three
     * strings to the console: "Phase0", "Integrator0", and "Phase1".
     */
    public static void main(String[] argv) {
        System.out.println(NameMaker.makeName(Phase.class));
        System.out.println(NameMaker.makeName(Integrator.class));
        System.out.println(NameMaker.makeName(Phase.class));
    }
}
