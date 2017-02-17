package etomica.config;

import etomica.api.IAtomList;
import etomica.api.IVector;

public class ConformationGeneric implements IConformation {

	protected final IVector[] coords;
	
	public ConformationGeneric(IVector[] coords) {
		this.coords = coords;
	}
	
	public void initializePositions(IAtomList atomList) {
		for (int i=0; i<atomList.getAtomCount(); i++) {
			atomList.getAtom(i).getPosition().E(coords[i]);
		}
	}

}
