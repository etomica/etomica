/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Grows configurations of chains with the pair probability distribution that
 * is flat for some distance and then has a power-law decay tail.
 */
public class MCMoveClusterAtomChainHSTail extends MCMoveBoxStep {

    protected final IRandom random;
    protected final double sigma;
    protected final Vector dr;
    protected final double pow;
    protected final double totPCore;

    public MCMoveClusterAtomChainHSTail(IRandom random, Box box, double sigma, double pow) {
        super();
        setBox(box);
        this.random = random;
        this.sigma = sigma;
        dr = box.getSpace().makeVector();
        this.pow = pow;
        
        totPCore = (pow-3) / pow;
    }

    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.size();

        for (int i=1; i<n; i++) {
            Vector pos = leafAtoms.get(i).getPosition();
            double rand = random.nextDouble();
            if (rand < totPCore) {
                pos.setRandomInSphere(random);
                pos.TE(sigma);
            }
            else {
                double r = Math.pow(1- (rand - totPCore)/(1-totPCore), -1.0/(pow-3))*sigma;
                pos.setRandomSphere(random);
                pos.TE(r);
            }
            pos.PE(leafAtoms.get(i-1).getPosition());
        }

		((BoxCluster)box).trialNotify();
		return true;
	}

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }
}
