package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomTreeNodeFactory;
import etomica.atom.AtomTreeNodeGroup;

/**
 * 
 */
public class AtomTreeNodeWater extends AtomTreeNodeGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeWater(Atom atom) {
		super(atom);
	}
	
	public Atom firstLeafAtom() {return O;}
	public Atom lastLeafAtom() {return H2;}

	public Atom H1, H2, O;
	
	public static final AtomTreeNodeFactory FACTORY = new AtomTreeNodeFactory() {
		public etomica.atom.AtomTreeNode makeNode(Atom atom) {
			return new AtomTreeNodeWater(atom);
		}
	};
	
	
}
