/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

public interface BetaSource {

    /**
     * Returns k such that beta[j] == 0 for all j > k. Largest index for a non-zero beta. 
     */
    public int maxIndex();
    
    public double beta(int k);
    
    /**
     * Returns k * beta(k), which in some cases has better precision than beta(k) itself.
     */
    public double kbeta(int k);

}
