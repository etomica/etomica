/*
 * History
 * Created on Nov 17, 2004 by kofke
 */
package etomica.data;

import etomica.DataTranslator;
import etomica.Space;
import etomica.space.Vector;

/**
 * Converts between an array of double and a vector object.
 * All methods assume that array is of correct dimension for
 * the vector.
 */
public class DataTranslatorVector implements DataTranslator {

    /**
     * Creates a translator that converts between a double and
     * a vector defined for the given space.
     */
    public DataTranslatorVector(Space space) {
        vector = space.makeVector();
        array = new double[space.D()];
    }

    
    public Object fromArray(double[] x) {
        vector.E(x);
        return vector;
    }

    public double[] toArray(Object obj) {
        ((Vector)obj).assignTo(array);
        return array;
    }

    private final Vector vector;
    private final double[] array;
}
