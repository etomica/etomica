/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Grows configurations of chains with the pair probability distribution that
 * is flat for some distance and then has a power-law decay tail.
 */
public class MCMoveClusterAtomChainHSTail extends MCMoveAtom {

    public MCMoveClusterAtomChainHSTail(IRandom random, ISpace _space, double sigma, double pow) {
        super(random, null, _space);
        this.sigma = sigma;
        dr = (IVectorRandom)space.makeVector();
        this.pow = pow;
        
        totPCore = (pow-3) / pow;
    }

    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.getAtomCount();

        for (int i=1; i<n; i++) {
            IVectorRandom pos = (IVectorRandom)leafAtoms.getAtom(i).getPosition();
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
	
    public double getA() {
        return 1;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected final double sigma;
    protected final IVectorRandom dr;
    protected final double pow;
    protected final double totPCore;
}
