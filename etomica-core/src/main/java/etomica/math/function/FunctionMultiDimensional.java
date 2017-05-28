/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.function;


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
