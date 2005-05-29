package etomica.models.water;

import etomica.Atom;
import etomica.AtomTreeNodeFactory;
import etomica.AtomTreeNodeGroup;

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
		public etomica.AtomTreeNode makeNode(Atom atom) {
			return new AtomTreeNodeWater(atom);
		}
	};
	
	
}
