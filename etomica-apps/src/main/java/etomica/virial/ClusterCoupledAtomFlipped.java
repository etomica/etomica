/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.IMolecule;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;

public class ClusterCoupledAtomFlipped  extends ClusterCoupledFlipped {

    protected final Vector axis;
    
    public ClusterCoupledAtomFlipped(ClusterAbstract cluster, Space space) {
        this(cluster,space, 0);
    }

    public ClusterCoupledAtomFlipped(ClusterAbstract cluster, Space space, double minFlipDistance) {
    	super(cluster,space, minFlipDistance);
    	axis = space.makeVector();
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledAtomFlipped(wrappedCluster.makeCopy(), space, minFlipDistance);
    }
   
    // flip the atom, not the molecule 
    protected void flip(IMolecule flippedMolecule) {
    	IAtomList childAtoms = flippedMolecule.getChildList();
    	if (childAtoms.size() > 1){
            throw new RuntimeException("i can just handle single atom!");
    	}
    	IOrientation or = ((IAtomOriented)childAtoms.get(0)).getOrientation();
    	if (or instanceof OrientationFull3D) {
    	    ((OrientationFull3D) or).getSecondaryDirection().normalize();
    	    axis.E(((OrientationFull3D) or).getSecondaryDirection());
    	    ((OrientationFull3D)or).rotateBy(Math.PI, axis);
    	    if (false && (Math.abs(axis.squared()-1) > 1e-10 || Math.abs(or.getDirection().squared() - 1) > 1e-10 || Math.abs(((OrientationFull3D)or).getSecondaryDirection().squared() - 1) > 1e-10)) {
    	        throw new RuntimeException("oops "+axis.squared()+" "+or.getDirection().squared()+" "+((OrientationFull3D)or).getSecondaryDirection().squared());
    	    }
    	}
    	else {
            Vector v = or.getDirection();
            v.TE(-1.0);
    	}
    }

}
