/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.atom.IAtomList;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;

/**
 * Sums the force and torque on each iterated atom or molecule and adds it to
 * the agent associated with the atom.
 */
public class PotentialCalculationForceSumWallForce extends PotentialCalculationForceSum implements IPotentialCalculationWallForce {

    public PotentialCalculationForceSumWallForce(IPotential potentialTether) {
        this.potentialTether = potentialTether;
    }
    
    /**
     * Re-zeros the force vectors.
     *
     */
    public void reset(){
        wallForce = 0;
        super.reset();
    }
    
    public double getWallForce() {
        return 0.5*wallForce;
    }
    
    /**
     * Adds forces and torques due to given potential acting on the atoms produced by the iterator.
     * Implemented for 1-, 2- and N-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        int nBody = potential.nBody();
        if (potential instanceof PotentialSoft) {
            PotentialSoft potentialSoft = (PotentialSoft)potential;
            Vector[] gradient = potentialSoft.gradient(atoms);
            switch(nBody) {
                case 1:
                    if (potential == potentialTether) {
                        if (atoms.getAtom(0).getPosition().getX(0) > 0) {
                            wallForce += gradient[0].getX(0);
                        }
                        else {
                            wallForce -= gradient[0].getX(0);
                        }
                    }
                    integratorAgentManager.getAgent(atoms.getAtom(0)).force().ME(gradient[0]);
                    break;
                case 2:
                    integratorAgentManager.getAgent(atoms.getAtom(0)).force().ME(gradient[0]);
                    integratorAgentManager.getAgent(atoms.getAtom(1)).force().ME(gradient[1]);
                    break;
                default:
                    //XXX atoms.count might not equal f.length.  The potential might size its 
                    //array of vectors to be large enough for one AtomSet and then not resize it
                    //back down for another AtomSet with fewer atoms.
                    for (int i=0; i<atoms.getAtomCount(); i++) {
                        integratorAgentManager.getAgent(atoms.getAtom(i)).force().ME(gradient[i]);
                    }
            }
        }
    }

    private static final long serialVersionUID = 1L;
	protected final IPotential potentialTether;
	protected double wallForce;
}
