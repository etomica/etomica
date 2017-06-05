package etomica.config;

import etomica.atom.IAtomList;
import etomica.space.Vector;

public class ConformationGeneric implements IConformation {

	protected final Vector[] coords;
	
	public ConformationGeneric(Vector[] coords) {
		this.coords = coords;
	}
	
	public void initializePositions(IAtomList atomList) {
		for (int i=0; i<atomList.getAtomCount(); i++) {
			atomList.getAtom(i).getPosition().E(coords[i]);
		}
	}

}
