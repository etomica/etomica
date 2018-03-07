/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.integrator.mcmove.MCMoveRotate;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveRotateNoSmer extends MCMoveRotate {

	public MCMoveRotateNoSmer(PotentialMaster potentialMaster, IRandom random,
                              Space _space) {
		super(potentialMaster, random, _space);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
	}

    public double getChi(double temperature) {
        if (associationManager.getAssociatedAtoms(atom).size() > 1) {
        	return 0;
        }
        return super.getChi(temperature);
    }
	protected AssociationManager associationManager;

}
