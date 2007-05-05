package etomica.models.water;

import etomica.atom.AtomGroup;
import etomica.atom.AtomType;
import etomica.atom.IAtomPositioned;

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
	
    private static final long serialVersionUID = 1L;
	public IAtomPositioned H1, H2, O;
}
