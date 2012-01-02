package etomica.virial.GUI.components;

import etomica.potential.P3BondAngle;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;

public abstract class ACollectionBondedPotential {
	
	
	public abstract void  addBondedPotentialSets(CollectionPotentialAtomicLike pObject, ISpace space, int speciesIndex);

	

}
