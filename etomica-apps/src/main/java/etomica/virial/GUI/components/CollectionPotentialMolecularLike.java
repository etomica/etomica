/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import etomica.potential.IPotentialMolecular;

public class CollectionPotentialMolecularLike implements ICollectionPotential{
	
	private int speciesIndex;
	//Array Object to declare Molecular Potential - Pairwise Interactions
	private IPotentialMolecular potentialMolecular;

	public  CollectionPotentialMolecularLike(int index){
		speciesIndex = index;
	}
	
	public void setPotentialMolecularNonBondedLike(IPotentialMolecular potentialMolecular) {
		// TODO Auto-generated method stub
		this.potentialMolecular = potentialMolecular;
	}

	public IPotentialMolecular getPotentialMolecularNonBondedLike() {
		// TODO Auto-generated method stub
		return this.potentialMolecular;
	}

	public int getSpeciesIndex() {
		return speciesIndex;
	}

	public void setSpeciesIndex(int speciesIndex) {
		this.speciesIndex = speciesIndex;
	}
	
	
}
