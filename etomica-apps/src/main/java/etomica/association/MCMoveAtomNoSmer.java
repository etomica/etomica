/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveAtomNoSmer extends MCMoveAtom {

	public MCMoveAtomNoSmer(IRandom random, PotentialMaster potentialMaster,
			Space _space) {
		super(random, potentialMaster, _space);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
	}

    public double getChi(double temperature) {
        if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 1) {
        	return 0;
        }
        return super.getChi(temperature);
    }
	protected AssociationManager associationManager;
}
