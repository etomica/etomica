package etomica.virial.GUI.components;

import etomica.potential.P3BondAngle;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;

public abstract class BondedInteraction {
	
	
	public void AddBondedPotentialSets(PotentialObject pObject, ISpace space, int speciesIndex){}

}
