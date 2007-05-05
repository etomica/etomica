package etomica.models.propane;

import etomica.atom.AtomGroup;
import etomica.atom.AtomType;
import etomica.atom.IAtomPositioned;

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
	
    private static final long serialVersionUID = 1L;
	public IAtomPositioned UA1, UA2, UA3;
}
