/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.OrientationFull3D;

public class ClusterCoupledAtomFlipped  extends ClusterCoupledFlipped {

    protected final IVectorMutable axis;
    
    public ClusterCoupledAtomFlipped(ClusterAbstract cluster, ISpace space) {
        this(cluster,space, 0);
    }

    public ClusterCoupledAtomFlipped(ClusterAbstract cluster, ISpace space, double minFlipDistance) {
    	super(cluster,space, minFlipDistance);
    	axis = space.makeVector();
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledAtomFlipped(wrappedCluster.makeCopy(), space, minFlipDistance);
    }
   
    // flip the atom, not the molecule 
    protected void flip(IMolecule flippedMolecule) {
    	IAtomList childAtoms = flippedMolecule.getChildList();
    	if (childAtoms.getAtomCount() > 1){
            throw new RuntimeException("i can just handle single atom!");
    	}
    	IOrientation or = ((IAtomOriented)childAtoms.getAtom(0)).getOrientation();
    	if (or instanceof OrientationFull3D) {
    	    ((IVectorMutable)((OrientationFull3D) or).getSecondaryDirection()).normalize();
    	    axis.E(((OrientationFull3D) or).getSecondaryDirection());
    	    ((OrientationFull3D)or).rotateBy(Math.PI, axis);
    	    if (false && (Math.abs(axis.squared()-1) > 1e-10 || Math.abs(or.getDirection().squared() - 1) > 1e-10 || Math.abs(((OrientationFull3D)or).getSecondaryDirection().squared() - 1) > 1e-10)) {
    	        throw new RuntimeException("oops "+axis.squared()+" "+or.getDirection().squared()+" "+((OrientationFull3D)or).getSecondaryDirection().squared());
    	    }
    	}
    	else {
            IVectorMutable v = (IVectorMutable)or.getDirection();
            v.TE(-1.0);
    	}
    }

}
