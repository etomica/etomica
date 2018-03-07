/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

import etomica.association.BiasVolume2SiteAceticAcid;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialGroup;

public class P2AceticAcidTwoSite extends PotentialGroup {
	double minE;
	protected boolean isAssociation;
	protected BiasVolume2SiteAceticAcid bv;
	protected int bondType;
    
	public P2AceticAcidTwoSite(double minE, boolean isAssociation, int bondType){
		super(2);
		this.isAssociation = isAssociation;
		this.bondType = bondType;
		this.minE = minE;
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
		boolean geometryAB = true;
		boolean geometryBA = true;
		IMolecule molecule1 = molecules.get(0);
		IMolecule molecule2 = molecules.get(1);

		if (bv != null){
			geometryAB = bv.isAssociated(molecule1, molecule2);
			geometryBA = bv.isAssociated(molecule2, molecule1);
	        if (bondType == 0 && !geometryAB){
	        	return 0.0;
	        }
	        if (bondType == 1 && !geometryBA){
	        	return 0.0;
	        }
		}

		double energy = super.energy(molecules); 

		if ((energy > -minE) == (bondType != 3)){//no association
			return 0.0;
    	}
		if (geometryAB && geometryBA){
        	energy = 0.5 * energy;
        }

		return energy;
	}
}
