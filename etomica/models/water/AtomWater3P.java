package etomica.models.water;

import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

/**
 * Tree node for 3-point water molecule.
 */
public class AtomWater3P extends AtomGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomWater3P(AtomType waterType) {
		super(waterType);
	}
	
	public AtomLeaf H1, H2, O;
}
