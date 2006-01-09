package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeFactory;
import etomica.atom.AtomTreeNodeGroup;

/**
 * Tree node for 3-point water molecule.
 */
public class AtomTreeNodeWater3P extends AtomTreeNodeGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeWater3P(Atom atom) {
		super(atom);
	}
	
	public AtomLeaf firstLeafAtom() {return O;}
	public AtomLeaf lastLeafAtom() {return H2;}

	public AtomLeaf H1, H2, O;
	
	public static final AtomTreeNodeFactory3P FACTORY = new AtomTreeNodeFactory3P();
    
    public static class AtomTreeNodeFactory3P implements AtomTreeNodeFactory, java.io.Serializable {
		public etomica.atom.AtomTreeNode makeNode(Atom atom) {
			return new AtomTreeNodeWater3P(atom);
		}
	};
	
	
}
