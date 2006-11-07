package etomica.models.propane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeFactory;
import etomica.atom.AtomTreeNodeGroup;

/**
 * Tree node for 3-point water molecule.
 */
public class AtomTreeNodeUAPropane extends AtomTreeNodeGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeUAPropane(Atom atom) {
		super(atom);
	}
	
	public AtomLeaf firstLeafAtom() {return UA1;}
	public AtomLeaf nextLeafAtom() {return UA2;}
	public AtomLeaf lastLeafAtom() {return UA3;}

	public AtomLeaf UA1, UA2, UA3;
	
	public static final AtomTreeNodeFactoryUAPropane FACTORY = new AtomTreeNodeFactoryUAPropane();
    
    public static class AtomTreeNodeFactoryUAPropane implements AtomTreeNodeFactory, java.io.Serializable {
		public etomica.atom.AtomTreeNode makeNode(Atom atom) {
			return new AtomTreeNodeUAPropane(atom);
		}
	};
	
	
}
