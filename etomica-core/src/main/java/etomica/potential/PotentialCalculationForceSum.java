/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.storage.VectorStorage;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum implements PotentialCalculation {
        
    protected final VectorStorage forces;

    public PotentialCalculationForceSum(VectorStorage forces) {
        this.forces = forces;
    }

    /**
     * Re-zeros the force vectors.
     *
     */
    public void reset(){
        for (int i = 0; i < forces.size(); i++) {
            forces.get(i).E(0);
        }
    }

    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        PotentialSoft potentialSoft = (PotentialSoft)potential;
        int nBody = potential.nBody();
        Vector[] f = potentialSoft.gradient(atoms);
        if (f==null) return;
        switch(nBody) {
        	case 0:
        		for (int i = 0; i< forces.size(); i++) {
        			forces.get(i).ME(f[i]);
        		}
        		break;

            case 1:
                forces.get(atoms.get(0).getLeafIndex()).ME(f[0]);
                break;
            case 2:
                forces.get(atoms.get(0).getLeafIndex()).ME(f[0]);
                forces.get(atoms.get(1).getLeafIndex()).ME(f[1]);
                break;
            default:
                //TODO atoms.count might not equal f.length.  The potential might size its
                //array of vectors to be large enough for one IAtomSet and then not resize it
                //back down for another IAtomSet with fewer atoms.
                for (int i = 0; i<atoms.size(); i++) {
                    forces.get(atoms.get(i).getLeafIndex()).ME(f[i]);
                }
		}
	}
}
