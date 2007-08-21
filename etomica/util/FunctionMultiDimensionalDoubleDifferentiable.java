package etomica.util;


/**
 * Interface for multi-dimensional 2-D differentiable function
 * 
 * @author Tai Boon Tan
 *
 */
public interface FunctionMultiDimensionalDoubleDifferentiable extends FunctionMultiDimensionalDifferentiable {
	
	public double[][] d2fdx2(double[] variables);
}
