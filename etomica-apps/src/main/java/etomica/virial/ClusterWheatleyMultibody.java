/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams, but adds in non-additive contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibody extends ClusterWheatleySoft {

    protected final MayerFunctionNonAdditive fNonAdditive;
    protected final MayerFunctionNonAdditive[] fMulti;
    protected final int[] moleculeIndices;
    protected final double[] r2;
    protected final MoleculeArrayList molecules;
    protected boolean doMulti;
    protected double rCut2;
    protected ClusterWheatleyMultibodyBD clusterMultiBD;
    protected double multiTol;
    protected final double[] fQmulti;

    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fMulti array of non-additive Mayer functions.  fMulti[3] is the
     *          3-body Mayer function (exp(-beta*deltaU3)-1), fMulti[4] is the
     *          4-body Mayer function, etc.  fMulti null entries will be
     *          ignored and the array need to not be of size equal to nPoints.
     *          If only 3-body Mayer function is available, then fMulti can be
     *          of length 4 (0,1,2,3).
     */
    public ClusterWheatleyMultibody(int nPoints, MayerFunction f, MayerFunctionNonAdditive[] fMulti) {
        this(nPoints, f, null, fMulti, 1e-12);
    }

    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fNonAdditive Mayer function that returns non-additive value for
     *          any number of molecules.
     */
    public ClusterWheatleyMultibody(int nPoints, MayerFunction f, MayerFunctionNonAdditive fNonAdditive) {
        this(nPoints, f, fNonAdditive, new MayerFunctionNonAdditive[0], 1e-12);
    }

    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fNonAdditive Mayer function that returns non-additive value for
     *          any number of molecules.
     * @param fMulti array of non-additive Mayer functions.  fMulti[3] is the
     *          3-body Mayer function (exp(-beta*deltaU3)-1), fMulti[4] is the
     *          4-body Mayer function, etc.  fMulti null entries will be
     *          ignored and the array need to not be of size equal to nPoints.
     *          If only 3-body Mayer function is available, then fMulti can be
     *          of length 4 (0,1,2,3).
     * @param tol if the magnitude of the computed cluster value is less than
     *          the value will be recomputed using BigDecimal.  Use tol=0 to
     *          prevent BigDecimal computations.
     */
    public ClusterWheatleyMultibody(int nPoints, MayerFunction f, MayerFunctionNonAdditive fNonAdditive, MayerFunctionNonAdditive[] fMulti, double tol) {
        super(nPoints, f, 0);
        this.fNonAdditive = fNonAdditive;
        this.fMulti = fMulti;
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // pairwise shouldn't try to use BD
        clusterBD = null;
        setTolerance(tol);
        molecules = new MoleculeArrayList(nPoints);
        rCut2 = Double.POSITIVE_INFINITY;
        fQmulti = new double[1<<n];
        fQmulti[0] = fQmulti[1] = fQmulti[2] = 1;
    }

    public void setTolerance(double newTol) {
        clusterMultiBD = new ClusterWheatleyMultibodyBD(n, f, fMulti, -3*(int)Math.log10(newTol));
        clusterMultiBD.setDoCaching(false);
        clusterMultiBD.setPrecisionLimit(300);
        multiTol = newTol;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleyMultibody c = new ClusterWheatleyMultibody(n, f, fNonAdditive, fMulti, multiTol);
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
        return c;
    }

    public void setTemperature(double newT) {
        super.setTemperature(newT);
        if (clusterMultiBD != null) {
            clusterMultiBD.setTemperature(newT);
        }
    }

    public void setRCut(double newRCut) {
        rCut2 = newRCut * newRCut;
    }

    public void calcValue(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        double rMax = 0;
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                if (cPairs.getr2(i,j) > rCut2) {
                    value = 0;
                    return;
                }
                if (cPairs.getr2(i,j) > rMax) rMax = cPairs.getr2(i,j);
            }
        }
        // do (multi+pair) - pair here so that we avoid recomputing f bonds
        doMulti = false;
        super.calcValue(box);
        double pairValue = value;
        doMulti = true;
        super.calcValue(box);
        value -= pairValue;
        // we have our own BD cluster for multibody
        // if an individual integrand (pair/total) is small, it won't trigger a BD calculation
        // BD only gets triggered here
        double bfac = (1.0-n)/SpecialFunctions.factorial(n);
        if (Math.abs(value) > 0 && Math.abs(value/bfac) < multiTol) {
            if (clusterMultiBD != null) {
                value = clusterMultiBD.value(box);
            }
            else {
                value = 0;
            }
        }
    }

    protected void calcFullFQ(BoxCluster box) {
        super.calcFullFQ(box);
        if (!doMulti) return;
        for (int i=3; i<fMulti.length; i++) {
            if (fMulti[i]!=null) fMulti[i].setBox(box);
        }
        if (fNonAdditive != null) {
            fNonAdditive.setBox(box);
        }
        int nf = 1<<n;
        IMoleculeList boxMolecules = box.getMoleculeList();
        // FQ[i] now contains the exp(-bU2) where U2 is the pair-wise energy for set i.
        // we need to go around and add the non-additive energy for each set.

        for (int i=3; i<nf; i++) {
            fQmulti[i] = 1;
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            if (fQ[i] == 0) continue;

            // we want to loop over subsets of i with at least 3 points
            // we could try to be clever and use the moleculeIndices constructed below to loop
            // more efficiently here.
            int iLowBit = (i & -i);//next lowest bit
            for (int isub=iLowBit; isub<i; isub+=iLowBit) {//sum over partitions of i
                while ((isub & ~i) != 0) {
                    // loop until isub is an actual subset of i
                    isub += iLowBit;
                }
                fQ[i] *= fQmulti[isub];
                if (fQ[i]==0) break;
            }
            if (fQ[i]==0) continue;
            

            int l = 0;
            molecules.clear();
            for (int a=0; a<n; a++) {
                if ((i & (1<<a)) != 0) {
                    moleculeIndices[l] = a;
                    molecules.add(boxMolecules.getMolecule(a));
                    l++;
                }
            }
            if ((fMulti.length <= l || fMulti[l] == null) && fNonAdditive == null) {
                fQmulti[i] = 1;
                continue;
            }
            int ll = 0;
            for (int a=0; a<l-1; a++) {
                for (int b=a+1; b<l; b++) {
                    r2[ll] = box.getCPairSet().getr2(moleculeIndices[a],moleculeIndices[b]);
                    ll++;
                }
            }
            if (fMulti.length > l && fMulti[l] != null) {
                fQmulti[i] = fMulti[l].f(molecules, l, moleculeIndices, r2, beta)+1;
            }
            fQ[i] *= fQmulti[i];
            if (fNonAdditive != null) {
                // we don't want to include this in fQmulti because we would just include it again
                // for larger sets
                fQ[i] *= fNonAdditive.f(molecules, l, moleculeIndices, r2, beta)+1;
            }
        }
    }
}
