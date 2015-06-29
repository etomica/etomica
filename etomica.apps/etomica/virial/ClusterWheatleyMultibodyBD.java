/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.math.BigDecimal;

import etomica.api.IMoleculeList;
import etomica.atom.MoleculeArrayList;
import etomica.math.SpecialFunctions;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams, but adds in non-additive contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibodyBD extends ClusterWheatleySoftBD {

    protected final MayerFunctionNonAdditive fMulti;
    protected final int[] moleculeIndices;
    protected final double[] r2;
    protected final MoleculeArrayList molecules;
    protected boolean doMulti;
    protected double rCut2;

    public ClusterWheatleyMultibodyBD(int nPoints, MayerFunction f, MayerFunctionNonAdditive fMulti, int precision) {
        super(nPoints, f, precision);
        this.fMulti = fMulti;
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // XXX we don't have a clusterBD for multibody.
        // set it to null so the failure isn't silent
        molecules = new MoleculeArrayList(nPoints);
        rCut2 = Double.POSITIVE_INFINITY;
    }
    
    public ClusterAbstract makeCopy() {
        ClusterWheatleyMultibodyBD c = new ClusterWheatleyMultibodyBD(n, f, fMulti, mc.getPrecision());
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
        return c;
    }

    public void setRCut(double newRCut) {
        rCut2 = newRCut * newRCut;
    }

    public void calcValue(BoxCluster box) {
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
        int nf = 1<<n;
        BigDecimal fBPair = fB[nf-1];
        double pairValue = value;
        doMulti = true;
        super.calcValue(box);
//        System.out.println("BDpair "+pairValue+"  full "+value+"   NA "+(value-pairValue));
        value = (1-n)*fB[nf-1].subtract(fBPair, mc).doubleValue()/SpecialFunctions.factorial(n);
    }

    protected void calcFullFQ(BoxCluster box) {
        super.calcFullFQ(box);
        if (!doMulti) return;
        fMulti.setBox(box);
        int nf = 1<<n;
        IMoleculeList boxMolecules = box.getMoleculeList();
        // FQ[i] now contains the exp(-bU2) where U2 is the pair-wise energy for set i.
        // we need to go around and add the non-additive energy for each set.

        for (int i=3; i<nf; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            if (fQ[i].doubleValue() == 0) continue;
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
                    r2[ll] = box.getCPairSet().getr2(moleculeIndices[a],moleculeIndices[b]);
                    ll++;
                }
            }
            fQ[i] = fQ[i].multiply(new BigDecimal(fMulti.f(molecules, l, moleculeIndices, r2, beta),mc).add(BDONE));
        }
    }
}
