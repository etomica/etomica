package etomica.models.water;

import etomica.Atom;
import etomica.AtomTreeNodeGroup;
import etomica.AtomTreeNode;

/**
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AtomTreeNodeWater extends AtomTreeNodeGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeWater(Atom atom, AtomTreeNodeGroup parent) {
		super(atom, parent);
	}
	
	public Atom firstLeafAtom() {return O;}
	public Atom lastLeafAtom() {return H2;}
	public boolean isResizable() {return false;}

	public Atom H1, H2, O;
	
	public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
		public etomica.AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup parent) {
			return new AtomTreeNodeWater(atom, parent);
		}
	};
	
	
}
