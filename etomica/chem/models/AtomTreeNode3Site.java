package etomica.chem.models;

import etomica.Atom;
import etomica.AtomTreeNodeGroup;
import etomica.AtomTreeNode;

/**
 * A group node having three child atoms, which for convenience can be
 * referenced by the fields atom1, atom2, and atom3, respectively.
 */

//maybe should use an atom array?
public class AtomTreeNode3Site extends AtomTreeNodeGroup {

	/**
	 * Constructor simply invokes superclass constructor.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNode3Site(Atom atom, AtomTreeNodeGroup parent) {
		super(atom, parent);
	}
	
	public Atom atom1, atom2, atom3 ;
	
	public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
		public etomica.AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup parent) {
			return new AtomTreeNode3Site(atom, parent);
		}
	};
	
	/**
	 * Calls superclass method and assigns the atom to one of the atom
	 * references.
	 */
	public void addAtomNotify(Atom atom) {
		super.addAtomNotify(atom);
		if(atom1 == null) atom1 = atom;
		else if(atom2 == null) atom2 = atom;
		else if(atom3 == null) atom3 = atom;
		else throw new RuntimeException("Error in attempting to add too many children to node that expects only three");
	}
	
	

	/**
	 * @see etomica.AtomTreeNodeGroup#removeAtomNotify(etomica.Atom)
	 */
	public void removeAtomNotify(Atom atom) {
		super.removeAtomNotify(atom);
		if(atom1 == atom) atom1 = null;
		else if(atom2 == atom) atom2 = null;
		else if(atom3 == atom) atom3 = null;
		else throw new RuntimeException("Error in attempting to remove atom from node");
	}

}





