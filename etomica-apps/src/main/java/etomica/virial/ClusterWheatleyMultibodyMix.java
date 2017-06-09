/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IMoleculeList;
import etomica.atom.MoleculeArrayList;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams for a mixture, but adds in 3-body contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibodyMix extends ClusterWheatleySoftMix {

    protected final MayerFunctionNonAdditive[][][] mixFMulti3;
    protected final MayerFunctionNonAdditive[] fMap3;
    protected final int[] moleculeIndices;
    protected final double[] r2;
    protected final MoleculeArrayList molecules;
    protected boolean doMulti;
    protected boolean nonAdditiveOnly;
    protected double rCut2;

    /**
     * Constructs a cluster capable of handling mixtures of nonadditive
     * molecules.
     * 
     * @param nPoints  total # of molecules
     * @param nTypes   # of molecules of each type
     * @param f        2D array of pair Mayer functions
     *                 (0-0, 0-1, 0-2..., 0-n),
     *                 (1-0, 1-1, 1-2..., 1-n)
     *                 ...
     *                 (n-0, n-1, n-2..., n-n) 
     * @param fMulti3  3D array of 3-body Mayer functions
     */
    public ClusterWheatleyMultibodyMix(int nPoints, int[] nTypes, MayerFunction[][] f, MayerFunctionNonAdditive[][][] fMulti3, double tol, boolean nonAddOnly) {
        super(nPoints, nTypes, f, tol);
        nonAdditiveOnly = nonAddOnly;
        mixFMulti3 = fMulti3;
        fMap3 = new MayerFunctionNonAdditive[1<<nPoints];
        int iType = 0, jType = 0, kType = 0;
        int iSum = nTypes[0], jSum = 0, kSum = 0;
        // this code could handle larger sets of non-additive interactions if only 
        // we had a way to construct maps
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
                
                kSum = jSum;
                kType = jType;
                for (int k=j+1; k<nPoints; k++) {
                    while (k>=kSum) {
                        kType++;
                        kSum += nTypes[kType];
                    }
                    fMap3[(1<<i)|(1<<j)|(1<<k)] = fMulti3[iType][jType][kType];
                }
            }
        }
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // XXX we don't have a clusterBD for multibody.
        // set it to null so the failure isn't silent
        clusterBD = null;
        molecules = new MoleculeArrayList(nPoints);
        rCut2 = Double.POSITIVE_INFINITY;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleyMultibodyMix c = new ClusterWheatleyMultibodyMix(n, nTypes, mixF, mixFMulti3, tol, nonAdditiveOnly);
        c.setTemperature(1/beta);
        return c;
    }

    public void calcValue(BoxCluster box) {
        if (!nonAdditiveOnly) {
            doMulti = true;
            super.calcValue(box);
            return;
        }
        CoordinatePairSet cPairs = box.getCPairSet();
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                if (cPairs.getr2(i,j) > rCut2) {
                    value = 0;
                    return;
                }
            }
        }
        // do (multi+pair) - pair here so that we avoid recomputing f bonds
        doMulti = false;
        super.calcValue(box);
        double pairValue = value;
        doMulti = true;
        super.calcValue(box);
        value -= pairValue;
    }
    
    protected void calcFullFQ(BoxCluster box) {
        int nf = 1<<n;
        IMoleculeList boxMolecules = box.getMoleculeList();
        // generate all partitions and compute product of e-bonds for all pairs in partition
        CoordinatePairSet cpairs = box.getCPairSet();
        for (int i=3; i<nf; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            int jj = k & -k; // 2nd lowest bit
            if (k == jj) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            fQ[i] = fQ[k]; //initialize with previously-computed product of all pairs in partition, other than j
            if (fQ[i] == 0) continue;
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                fQ[i] *= fQ[l | j];
            }
            if (!doMulti || fMap3[i] == null || fQ[i] == 0) continue;

            // this code could easily handle larger sets (4, 5, etc)
            int kk = k&~jj; // strip jj bit from k
            int jjj = kk & -kk; // 3rd lowest bit
            if (kk == jjj) {
                // i is a 3-point set.  add in 3-body contribution
                int l = 0;
                molecules.clear();
                for (int a=0; a<n; a++) {
                    if ((i & (1<<a)) != 0) {
                        moleculeIndices[l] = a;
                        molecules.add(boxMolecules.getMolecule(a));
                        l++;
                    }
                }
                int ll = 0;
                for (int a=0; a<l-1; a++) {
                    for (int b=a+1; b<l; b++) {
                        r2[ll] = cpairs.getr2(moleculeIndices[a],moleculeIndices[b]);
                        ll++;
                    }
                }
                fMap3[i].setBox(box);
                fQ[i] *= fMap3[i].f(molecules, 3, moleculeIndices, r2, beta)+1;
            }
        }
    }
}
