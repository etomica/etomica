/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.overlap;

/**
 * Class that can return values of alpha.
 * 
 * @author Andrew Schultz
 */
public interface AlphaSource {
    
    /**
     * Returns the number of possible alpha values.
     */
    public int getNumAlpha();

    /**
     * Returns the ith alpha value (i=0..(n-1))
     */
    public double getAlpha(int i);
}