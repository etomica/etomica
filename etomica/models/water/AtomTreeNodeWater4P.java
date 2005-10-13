package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomTreeNodeFactory;
import etomica.atom.AtomTreeNodeGroup;

/**
 * Tree node for a 4-point water molecule.
 */
public class AtomTreeNodeWater4P extends AtomTreeNodeGroup {

	/**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomTreeNodeWater4P(Atom atom) {
		super(atom);
	}
	
	public Atom firstLeafAtom() {return O;}
	
	public Atom lastLeafAtom() {return H2;}

	public Atom H1, H2, O, M;
	
    public static final AtomTreeNodeFactory4P FACTORY = new AtomTreeNodeFactory4P();
    
    public static class AtomTreeNodeFactory4P implements AtomTreeNodeFactory, java.io.Serializable {
        public etomica.atom.AtomTreeNode makeNode(Atom atom) {
            return new AtomTreeNodeWater4P(atom);
        }
    };
}
