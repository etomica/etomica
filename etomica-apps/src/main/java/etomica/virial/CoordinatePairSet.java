/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


public interface CoordinatePairSet {

    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j);

    /**
     * Informs the CoordinatePairSet that the configuration has changed and that it
     * has a new ID
     */
    public void reset(long cPairID);

    public long getID();

}