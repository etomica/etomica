/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IPotentialAtomic;
import etomica.potential.Potential;
import etomica.space.Space;

public class PotentialThreaded extends Potential {

	final protected IPotentialAtomic[] potential;
	
	public PotentialThreaded(Space space, IPotentialAtomic[] potential) {
		super(potential[0].nBody(), space);
		this.potential = potential;
	
	}

	public double energy(IAtomList atoms) {
		//Only the energy from one thread (a partition of atoms)
		return potential[0].energy(atoms);
	}

	public double getRange() {
		
		return potential[0].getRange();
	}

	public void setBox(Box box) {
		
		for(int i=0; i<potential.length; i++){
			potential[i].setBox(box);
		}

	}
	
	public IPotentialAtomic[] getPotentials(){
		return potential;
	}

}
