package etomica.models.water;

import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

/**
 * Tree node for a 4-point water molecule.
 */
public class AtomWater4P extends AtomGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomWater4P(AtomType waterType) {
		super(waterType);
	}
	
	public AtomLeaf H1, H2, O, M;
}
