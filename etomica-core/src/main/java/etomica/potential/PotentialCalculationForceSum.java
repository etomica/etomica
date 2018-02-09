/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorBox;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum implements PotentialCalculation {
        
    protected AtomLeafAgentManager<Vector> integratorAgentManager;

    public void setAgentManager(AtomLeafAgentManager<Vector> agentManager) {
        integratorAgentManager = agentManager;
    }

    /**
     * Re-zeros the force vectors.
     *
     */
    public void reset(){
        this.integratorAgentManager.getAgents().values().forEach((agent) -> agent.E(0));
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
        		IAtomList boxAtoms = integratorAgentManager.getBox().getLeafList();
        		for (int i=0; i<boxAtoms.getAtomCount(); i++) {
        			integratorAgentManager.getAgent(boxAtoms.getAtom(i)).ME(f[i]);
        		}
        		break;

            case 1:
                integratorAgentManager.getAgent(atoms.getAtom(0)).ME(f[0]);
                break;
            case 2:
                integratorAgentManager.getAgent(atoms.getAtom(0)).ME(f[0]);
                integratorAgentManager.getAgent(atoms.getAtom(1)).ME(f[1]);
                break;
            default:
                //TODO atoms.count might not equal f.length.  The potential might size its
                //array of vectors to be large enough for one IAtomSet and then not resize it
                //back down for another IAtomSet with fewer atoms.
                for (int i=0; i<atoms.getAtomCount(); i++) {
                    integratorAgentManager.getAgent(atoms.getAtom(i)).ME(f[i]);
                }
		}
	}
}
