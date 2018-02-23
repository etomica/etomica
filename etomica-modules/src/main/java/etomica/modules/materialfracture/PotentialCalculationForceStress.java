/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;

import etomica.atom.IAtomList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialCalculationForceStress extends
        PotentialCalculationForcePressureSum {

    public PotentialCalculationForceStress(Space space) {
        super(space);
    }
    
    public void reset() {
        super.reset();
        load = 0;
    }

    public double getLoad() {
        return load;
    }

    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        PotentialSoft potentialSoft = (PotentialSoft)potential;
        int nBody = potential.nBody();
        Vector[] f = potentialSoft.gradient(atoms, pressureTensor);
        switch(nBody) {
            case 1:
                integratorAgentManager.getAgent(atoms.get(0)).ME(f[0]);
                if (potential instanceof P1Tension) {
                    load += Math.abs(f[0].getX(0));
                }
                break;
            case 2:
                integratorAgentManager.getAgent(atoms.get(0)).ME(f[0]);
                integratorAgentManager.getAgent(atoms.get(1)).ME(f[1]);
                break;
            default:
                //XXX atoms.count might not equal f.length.  The potential might size its 
                //array of vectors to be large enough for one IAtomSet and then not resize it
                //back down for another IAtomSet with fewer atoms.
                for (int i = 0; i<atoms.size(); i++) {
                    integratorAgentManager.getAgent(atoms.get(i)).ME(f[i]);
                }
        }
    }

    protected double load;
}
