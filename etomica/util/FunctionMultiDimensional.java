package etomica.util;


/**
 * Interface for multi-dimensional function.
 * 
 * @author Tai Boon Tan
 *
 */
public interface FunctionMultiDimensional {

	public double f(double[] variables);
    
    /**
     * The dimension of the space of independent variables.  The length of the
     * array passed to the f should be equal to this quantity. 
     */
    public int getDimension();

}
