/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;

public class MCMoveAtomNoSmer extends MCMoveAtom {

	public MCMoveAtomNoSmer(IRandom random, PotentialMaster potentialMaster,
			Space _space) {
		super(random, potentialMaster, _space);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
	}
	public double getA() {
        if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 1) {
        	return 0;
        } 
        return 1;
	}
	protected AssociationManager associationManager;
}
