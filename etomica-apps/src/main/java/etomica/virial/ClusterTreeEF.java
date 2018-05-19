/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Sep 24, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.virial.cluster.ClusterDiagramTree;

import java.util.Arrays;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterTreeEF extends ClusterTree {

    /**
     * @param bonds
     * @param eArray
     */
    public ClusterTreeEF(ClusterDiagramTree bonds, MayerFunction[] eArray) {
        super(bonds, Arrays.copyOf(eArray, eArray.length * 2));
        numF = f.length/2;
    }

    public ClusterAbstract makeCopy() {
        MayerFunction[] eArray = new MayerFunction[numF];
        System.arraycopy(f,0,eArray,0,numF);
        ClusterTreeEF copy = new ClusterTreeEF(bondsTree,eArray);
        copy.setTemperature(1/beta);
        return copy;
    }

    protected void updateF(CoordinatePairLeafSet cPairs, AtomPairSet aPairs) {
        int nPoints = pointCount();
        
        // recalculate all f values for all pairs
        int l=0;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<numF; k++) {
                    fValues[l][k+numF] = f[k].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                    fValues[l][k] = fValues[l][k+numF] - 1;
                }
                l++;
            }
        }
    }

    private final int numF;
}
