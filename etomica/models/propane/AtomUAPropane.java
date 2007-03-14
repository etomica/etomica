package etomica.models.propane;

import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

/**
 * Tree node for 3-point water molecule.
 */
public class AtomUAPropane extends AtomGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomUAPropane(AtomType moleculeType) {
		super(moleculeType);
	}
	
	public AtomLeaf UA1, UA2, UA3;
}
