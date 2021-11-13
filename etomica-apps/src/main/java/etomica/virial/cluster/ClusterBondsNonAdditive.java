/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

/**
 * Subclass of ClusterBonds capable of handling non-additive Mayer functions.
 * 
 * @author Andrew Schultz
 */
public class ClusterBondsNonAdditive extends ClusterBonds {

    /**
     * Construct a ClusterBondsNonAdditive of size n, using the given list of
     * pair and multibonds.  The pairBonds are indexed by i and j pairs, and
     * then by type of bond.  The multibonds are indexed by size of molecular
     * group, and then by groupID.
     */
    public ClusterBondsNonAdditive(int n, int[][][] pairBonds, int[][] multiBonds) {
        super(n, pairBonds, false);
        this.multiBonds = multiBonds;
    }

    public ClusterBondsNonAdditive(int[][] bondIndexArray, int[][] multiBonds) {
        super(bondIndexArray, false);
        this.multiBonds = multiBonds;
    }

    public double value(double[][][] fValues, double[][] fMultiValues) {
        // let the superclass handle any pairwise additive parts of the diagram
        double v = super.value(fValues);
        // now include the multibody parts
        for (int i=3; i<multiBonds.length; i++) {
            for (int j=0; j<multiBonds[i].length; j++) {
                v *= fMultiValues[i][multiBonds[i][j]];
            }
        }
        return v;
    }

    /**
     * Returns the list of multiBonds to multiply.  The list is indexed first
     * by the number of molecules in a group and then by groupID.
     */
    public int[][] getMultiBonds() {
        return multiBonds;
    }

    protected final int[][] multiBonds;
}
