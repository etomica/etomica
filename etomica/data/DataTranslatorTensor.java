/*
 * History
 * Created on Nov 17, 2004 by kofke
 */
package etomica.data;

import etomica.DataTranslator;
import etomica.Space;
import etomica.space.Tensor;

/**
 * Converts between an array of double and a tensor object.
 */
public class DataTranslatorTensor implements DataTranslator {

    /**
     * Creates a translator that converts between a double and
     * a tensor defined for the given space.
     */
    public DataTranslatorTensor(Space space) {
        tensor = space.makeTensor();
        array = new double[space.D * space.D];
    }

    
    public Object fromArray(double[] x) {
        tensor.E(x);
        return tensor;
    }

    public double[] toArray(Object obj) {
        ((Tensor)obj).assignTo(array);
        return array;
    }

    private final Tensor tensor;
    private final double[] array;
}
