package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomOriented;
import etomica.atom.IAtomPositionDefinition;
import etomica.space.ISpace;

public class ClusterCoupledAtomFlipped  extends ClusterCoupledFlipped {

    public ClusterCoupledAtomFlipped(ClusterAbstract cluster, ISpace space) {
    	super(cluster,space);

    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledAtomFlipped(wrappedCluster.makeCopy(), space);
    }
   
    // flip the atom, not the molecule 
    protected void flip(IMolecule flippedMolecule) {
    	IAtomList childAtoms = flippedMolecule.getChildList();
    	if (childAtoms.getAtomCount() > 1){
            throw new RuntimeException("i can just handle single atom!");
    	}
    	IAtomOriented atom = (IAtomOriented)childAtoms.getAtom(0);
        IVectorMutable v = (IVectorMutable) atom.getOrientation().getDirection();
        v.TE(-1.0);
    }

}
