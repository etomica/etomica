/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.function;

/**
 * 
 * Interface for multi-dimensional 1-D differentiable function
 * 
 * @author Tai Boon Tan
 *
 */
public interface FunctionMultiDimensionalDifferentiable extends FunctionMultiDimensional{
	
    /**
     * Returns the derivative indicated by the array d evaluated at the given point x.  The array d
     * should be the same length as x. Each entry of d specifies the order of the derivative
     * taken with respect to the corresponding variable in x.  For example, for a function
     * on a 2-dimensional space, d = {2,0}, indicates d2f/dx[0]^2, while d = {1,1}, indicates
     * d2f/dx[0]dx[1]; d = {2,1} indicates d3f/dx[0]^2dx[1]; etc.
     */
	public double df(int[] d, double[] x);
	
}
