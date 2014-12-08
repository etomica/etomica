/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class MCMoveClusterAtomHSChain extends MCMoveAtom {

    public MCMoveClusterAtomHSChain(IRandom random, ISpace _space, double sigma) {
        super(random, null, _space);
        this.sigma = sigma;
        dr = (IVectorRandom)space.makeVector();
    }
    
    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.getAtomCount();
        if (seq == null) {
            seq = new int[n];
        }
        for (int i=0; i<n; i++) {
            seq[i] = i;
        }
        for (int i=0; i<n; i++) {
            int j = i+random.nextInt(n-i);
            int k = seq[j];
            seq[j] = seq[i];
            seq[i] = k;
        }
        leafAtoms.getAtom(seq[0]).getPosition().E(0);

        for (int i=1; i<n; i++) {
            IVectorRandom pos = (IVectorRandom)leafAtoms.getAtom(seq[i]).getPosition();

            pos.setRandomInSphere(random);
            pos.TE(sigma);
            pos.PE(leafAtoms.getAtom(seq[i-1]).getPosition());
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
    protected int[] seq;
}
