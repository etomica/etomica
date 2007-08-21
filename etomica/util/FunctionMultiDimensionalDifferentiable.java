package etomica.util;

/**
 * 
 * Interface for multi-dimensional 1-D differentiable function
 * 
 * @author Tai Boon Tan
 *
 */
public interface FunctionMultiDimensionalDifferentiable extends FunctionMultiDimensional{
	
	public double[] dfdx (double[] variables);
	
}
