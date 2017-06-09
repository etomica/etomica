/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.association.BiasVolume2SiteAceticAcid;
import etomica.potential.PotentialGroup;

public class P2HardAssociationRefAceticAcidTwoSite extends PotentialGroup {
	double minE;
	protected BiasVolume2SiteAceticAcid bv;
    
	public P2HardAssociationRefAceticAcidTwoSite(){
		super(2);
	}
	
	public void setBiasVolume(BiasVolume2SiteAceticAcid bv){
		this.bv = bv;
	}
	
	public void setBox(Box box){
		super.setBox(box);
		if (bv != null){
			bv.setBox(box);
		}
	}

	public double energy(IMoleculeList molecules){
		boolean geometry = true;
		IMolecule molecule1 = molecules.getMolecule(0);
		IMolecule molecule2 = molecules.getMolecule(1);
		double energy = 0.0;
		
		geometry = bv.isAssociated(molecule1, molecule2);

		if (geometry){
			energy = Double.POSITIVE_INFINITY;
		}
		return energy;
	}
}
