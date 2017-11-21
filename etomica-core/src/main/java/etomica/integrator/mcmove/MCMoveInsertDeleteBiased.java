/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.random.IRandom;

import java.util.Arrays;

/**
 * Biased MCMove for insertion and deletion.  The bias can be set by setting
 * the natural log of the bias (lnBias) for every number of atoms.  The range
 * of atom numbers that will be sampled can also be set via fixedN and maxDN.
 * The range goes from (fixedN-maxDN) to (fixedN+maxDN).  Moves that would
 * leave that range are still proposed, but when accepted, the move is
 * rejected internally such that the system will remain in the range.
 *
 * @author Andrew Schultz
 */
public class MCMoveInsertDeleteBiased extends MCMoveInsertDelete {

    protected int minN, maxN;
    protected double[] lnbias;

    public MCMoveInsertDeleteBiased(PotentialMaster potentialMaster,
                                    IRandom random, Space _space, int minN, int maxN) {
        super(potentialMaster, random, _space);
        lnbias = new double[0];
        this.minN = minN;
        this.maxN = maxN;
    }
    
    public void setBox(Box box) {
        super.setBox(box);
    }
    
    public void setMinMax(int newMinN, int newMaxN) {
        minN = newMinN;
        maxN = newMaxN;
    }

    public double getLnBias(int n) {
        return lnbias[n];
    }

    public void setLnBias(int n, double nBias) {
        if (lnbias.length < n+1) {
            int oldSize = lnbias.length;
            lnbias = Arrays.copyOf(lnbias, n + 1);
            for (int i=oldSize; i<lnbias.length; i++) {
                if (i==0) lnbias[0] = 0;
                else lnbias[i] = lnbias[i-1] + mu;
            }
        }
        lnbias[n] = nBias;
    }

    public void setMu(double newMu) {
        super.setMu(newMu);
        if (lnbias != null && lnbias.length > 0) {
            lnbias[0] = 0;
            for (int i=1; i<lnbias.length; i++) {
                lnbias[i] = lnbias[i-1] + mu;
            }
        }
    }

    public double getChi(double temperature) {
        if (insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();
        }
        else {
            uNew = 0.0;
        }
        double b = uOld - uNew;

        int numMolecules = box.getNMolecules(species);
        double a = box.getBoundary().volume() / numMolecules;
        return Math.exp(getLnBiasDiff()) * a * Math.exp(b / temperature);
    }
    
    public double getLnBiasDiff() {
        int numAtoms = box.getLeafList().getAtomCount();
        if (lnbias != null && lnbias.length < numAtoms+1) {
            int oldSize = lnbias.length;
            lnbias = Arrays.copyOf(lnbias, numAtoms+1);
            for (int i=oldSize; i<lnbias.length; i++) {
                if (i==0) lnbias[0] = 0;
                //should be beta*mu here!!!
                else lnbias[i] = lnbias[i-1] + mu;
            }
        }
        double bias = lnbias == null ? 1 : lnbias[numAtoms]-lnbias[numAtoms-1];
        return (insert ? bias : -bias);
    }

    public void myAcceptNotify() {
        super.acceptNotify();
    }
    
    public void myRejectNotify() {
        super.rejectNotify();
    }
    
    public final void rejectNotify() {
        myRejectNotify();
    }
    
    public void acceptNotify() {
        int numAtoms = box.getLeafList().getAtomCount();
        if ((numAtoms <= minN && !insert) || (numAtoms > maxN && insert)) {
            myRejectNotify();
        }
        else {
            myAcceptNotify();
        }
    }
}
