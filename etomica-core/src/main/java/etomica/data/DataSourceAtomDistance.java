/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.atom.IAtom;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Length;

public class DataSourceAtomDistance extends DataSourceScalar {
	
    public DataSourceAtomDistance(Space space) {
		super("interatomic distance", Length.DIMENSION);
		
		this.space = space;
	
		vector = space.makeVector(); // to avoid making the vector each time getData() is called

	}


	public double getDataAsScalar() {

        vector.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        
		return Math.sqrt(vector.squared());
	}

	
	public void setAtoms(IAtom atom1, IAtom atom2) {
		this.atom1 = atom1;
		this.atom2 = atom2;
	}
	

    private static final long serialVersionUID = 1L;
	protected final Space space;
	protected final Vector vector;
	protected IAtom atom1;
	protected IAtom atom2;
}
