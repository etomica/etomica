/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Grows configurations of chains with the pair probability distribution that
 * is flat for some distance and then has a power-law decay tail.
 */
public class MCMoveClusterAtomChainHSTail extends MCMoveAtom {

    public MCMoveClusterAtomChainHSTail(IRandom random, Space _space, double sigma, double pow) {
        super(random, null, _space);
        this.sigma = sigma;
        dr = space.makeVector();
        this.pow = pow;
        
        totPCore = (pow-3) / pow;
    }

    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.getAtomCount();

        for (int i=1; i<n; i++) {
            Vector pos = leafAtoms.getAtom(i).getPosition();
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
            pos.PE(leafAtoms.getAtom(i-1).getPosition());
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

    protected final double sigma;
    protected final Vector dr;
    protected final double pow;
    protected final double totPCore;
}
