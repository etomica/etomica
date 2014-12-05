/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.space.ISpace;

public class MCMoveRotateNoSmer extends MCMoveRotate {

	public MCMoveRotateNoSmer(IPotentialMaster potentialMaster, IRandom random,
			ISpace _space) {
		super(potentialMaster, random, _space);
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
