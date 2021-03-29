/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.math.BigDecimal;

/**
 * Cluster class using Whealtey's recursion with BigDecimal to handle mixtures.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleySoftMixBD extends ClusterWheatleySoftBD {

    protected final MayerFunction[][] mixF;
    protected final int[] nTypes;
    protected final MayerFunction[][] fMap;
    
    public ClusterWheatleySoftMixBD(int nPoints, int[] nTypes, MayerFunction[][] f, int precision) {
        super(nPoints, null, precision);
        this.nTypes = nTypes;
        mixF = f;
        fMap = new MayerFunction[nPoints][nPoints];
        int iType = 0, jType = 0;
        int iSum = nTypes[0], jSum = 0;
        for (int i=0; i<nPoints; i++) {
            while (i>=iSum) {
                iType++;
                iSum += nTypes[iType];
            }
            jSum = iSum;
            jType = iType;
            for (int j=i+1; j<nPoints; j++) {
                while (j>=jSum) {
                    jType++;
                    jSum += nTypes[jType];
                }
                fMap[i][j] = f[iType][jType];
            }
        }
                
    }
    
    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftMixBD c = new ClusterWheatleySoftMixBD(n, nTypes, mixF, mc.getPrecision());
        c.setTemperature(1/beta);
        return c;
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        for (int i=0; i<mixF.length; i++) {
            for (int j=0; j<mixF[i].length; j++) {
                mixF[i][j].setBox(box);
            }
        }

        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double ff = fMap[i][j].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)] = new BigDecimal(ff).add(BDONE, mc);
            }
        }
    }

}
