package etomica.atom;

import etomica.Atom;

/**
 * A group node which holds also an array for convenient referencing of the
 * atoms. Appropriate for groups in which it is not expected that atoms will
 * be repeatedly added and removed.
 */

/*
 * History
 * 03/02/04 (DAK) new
 */
public class AtomTreeNodeGroupArray extends AtomTreeNodeGroup {

	/**
	 * Constructor simply invokes superclass constructor.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeGroupArray(Atom atom) {
		super(atom);
	}
	
	public Atom[] childAtomArray = new Atom[0];
	
	public static final AtomTreeNodeFactory FACTORY = new AtomTreeNodeFactory() {
		public etomica.atom.AtomTreeNode makeNode(Atom atom) {
			return new AtomTreeNodeGroupArray(atom);
		}
	};
	
	
	
	/**
	 * Calls superclass method and assigns the atom to one of the atom
	 * references.
	 */
	public void addAtomNotify(Atom atom) {
		super.addAtomNotify(atom);
		if(arrayIndex(atom) != -1) throw new RuntimeException("Adding child to node where it is already present as child");
		Atom[] newArray = new Atom[childAtomArray.length+1];
		System.arraycopy(childAtomArray, 0, newArray, 0, childAtomArray.length);
		newArray[childAtomArray.length] = atom;
		childAtomArray = newArray;
	}
	
	private int arrayIndex(Atom a) {
		for(int i=0; i<childAtomArray.length; i++) if(childAtomArray[i] == a) return i;
		return -1;
	}
	

	/**
	 * @see etomica.atom.AtomTreeNodeGroup#removeAtomNotify(etomica.Atom)
	 */
	public void removeAtomNotify(Atom atom) {
		super.removeAtomNotify(atom);
		int idx = arrayIndex(atom);
		if(idx == -1) throw new RuntimeException("Error in attempting to remove atom from node");
		Atom[] newArray = new Atom[childAtomArray.length-1];
		int k = 0;
		for(int i=0; i<childAtomArray.length; i++) {
			if(i == idx) continue;
			newArray[k++] = childAtomArray[i];
		}
		childAtomArray = newArray;
	}

	/**
	 * @see etomica.atom.AtomTreeNodeGroup#childAtomArray()
	 */
	public Atom[] childAtomArray() {
		return childAtomArray;
	}

}
