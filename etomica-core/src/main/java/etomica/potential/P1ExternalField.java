/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.box.Box;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

public class P1ExternalField extends PotentialMolecular {

	protected DipoleSource dipoleSource;
	protected Vector externalField;
	
	
	public P1ExternalField( Space space) {
		super(1,space);
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}
	
	public void setExternalField(Vector newExternalField){
		externalField = newExternalField;
	}
	

	public double energy(IMoleculeList molecules) {
		Vector dr = dipoleSource.getDipole(molecules.getMolecule(0));
		double energy = -externalField.dot(dr);
//		System.out.println("E = " + externalField + " dipole = " + dr );
		return energy;
	}

	public void setBox(Box box) {
		
	}
	
	public void setDipoleSource(DipoleSource newDipoleSource) {
		dipoleSource = newDipoleSource;
	}

}
