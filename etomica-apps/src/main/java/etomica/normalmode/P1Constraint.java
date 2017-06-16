/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.potential.Potential1;
import etomica.space.Space;

public class P1Constraint extends Potential1{

	public P1Constraint(Space space, double neighborRadius, Box box, CoordinateDefinition coordinateDefinition) {
	    super(space);

	    siteIndex = box.getLeafList().getAtomCount();
	    latticeSite = space.makeVectorArray(siteIndex);

	    radiusInner = neighborRadius*neighborRadius/4;

	    //Lattice Site Assignment
	    for(int i=0; i<siteIndex; i++){
	        latticeSite[i] = coordinateDefinition.getLatticePosition(box.getLeafList().getAtom(i));
	    }
	}
	
	public double getInnerRadius(){
	    return radiusInner;
	}
	
	@Override
	public double energy(IAtomList atoms) {

	    Vector posAtom = atoms.getAtom(0).getPosition();

	    int atomIndex = atoms.getAtom(0).getLeafIndex();
	    double d = posAtom.Mv1Squared(latticeSite[atomIndex]);

	    if (d < radiusInner){
	        return 0;
	    }

	    return Double.POSITIVE_INFINITY;
	}

    private static final long serialVersionUID = 1L;
    private Vector[] latticeSite;
    private int siteIndex;
    private double radiusInner;
}
