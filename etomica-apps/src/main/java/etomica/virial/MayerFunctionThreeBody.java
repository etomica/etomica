/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.virial.cluster.VirialDiagrams;

/**
 * This class acts as a non-additive Mayer function for molecules that
 * interact with only a 3-body potential.
 * 
 * This class handles bookkeeping and relies on a subclass to know how to
 * actually calculate the potential.
 *
 * @author Andrew Schultz
 */
public abstract class MayerFunctionThreeBody implements MayerFunctionNonAdditive {

    protected double[] lastValue;
    protected int totalMolecules;
    protected int[] myMoleculeIndices = new int[0];

    public double f(IMoleculeList molecules, double[] r2, double beta) {
        // just construct info to pass to the other method
        int nMolecules = molecules.getMoleculeCount();
        for (int i=0; i<nMolecules; i++) {
            myMoleculeIndices[i] = molecules.getMolecule(i).getIndex();
        }
        return f(molecules, nMolecules, myMoleculeIndices, r2, beta);
    }

    /*
     * We assume that this will be called for all triplets before it is called
     * for anything else.  We'll calculate the higher order energies by summing
     * the triplet contributions
     */
    public double f(IMoleculeList molecules, int nMolecules, int[] moleculeIndices, double[] r2, double beta) {
        if (nMolecules<3) throw new RuntimeException("need at least 3");
        if (nMolecules==3) {
            int id0 = moleculeIndices[0];
            int id1 = moleculeIndices[1];
            int id2 = moleculeIndices[2];
            int tripletID = VirialDiagrams.tripletId(id0, id1, id2, totalMolecules);
            double betaU = -beta*energy(molecules, r2);
            lastValue[tripletID] = betaU;
            return calcF(betaU);
        }
        double betaUSum = 0;
        // compute sum of triplet energies for all triplets composed of our molecules
        for (int i=0; i<nMolecules-2; i++) {
            int id0 = moleculeIndices[i];
            for (int j=i+1; j<nMolecules-1; j++) {
                int id1 = moleculeIndices[j];
                int tripletIDij = VirialDiagrams.tripletId(id0, id1, id1, totalMolecules);
                for (int k=j+1; k<nMolecules; k++) {
                    int id2 = moleculeIndices[k];
                    int tripletID = tripletIDij + (id2 - id1);
                    betaUSum += lastValue[tripletID];
                }
            }
        }
        return calcF(betaUSum);
    }

    /**
     * Returns exp(x)-1 computed directly or using a series expansion
     */
    public double calcF(double x) {
        if (Math.abs(x) < 0.01) {
            return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        return Math.exp(x)-1;
    }
    
    protected abstract double energy(IMoleculeList molecules, double[] r2);

    public void setBox(Box box) {
        if (lastValue != null) return;
        IMoleculeList molecules = box.getMoleculeList();
        totalMolecules = molecules.getMoleculeCount();
        if (totalMolecules<3) throw new RuntimeException("need at least 3");
        int nTriplets = VirialDiagrams.tripletId(totalMolecules-3, totalMolecules-2, totalMolecules-1, totalMolecules)+1;
        lastValue = new double[nTriplets];
        if (totalMolecules > myMoleculeIndices.length) {
            myMoleculeIndices = new int[totalMolecules];
        }
    }
}
