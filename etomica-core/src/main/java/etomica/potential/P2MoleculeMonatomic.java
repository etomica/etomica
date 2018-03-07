/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomPair;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;

public class P2MoleculeMonatomic implements IPotentialMolecular {
	public P2MoleculeMonatomic(IPotentialAtomic potential) {
		wrappedPotential = potential;
		leafAtoms = new AtomPair();
	}

	public double energy(IMoleculeList atoms) {
		leafAtoms.atom0 = atoms.get(0).getChildList().get(0);
		leafAtoms.atom1 = atoms.get(1).getChildList().get(0);
		return wrappedPotential.energy(leafAtoms); 
	}

	public double getRange() {
		return wrappedPotential.getRange();
	}

	public int nBody() {
		return 2;
	}

	public void setBox(Box box) {
		wrappedPotential.setBox(box); 

	}

	public IPotentialAtomic getWrappedPotential() { 
		return wrappedPotential; 
		} 
	public void setWrappedPotential(IPotentialAtomic newWrappedPotential) { 
		wrappedPotential = newWrappedPotential; 
		} 
	private static final long serialVersionUID = 1L; 
	protected final AtomPair leafAtoms;
	protected IPotentialAtomic wrappedPotential; 

}
