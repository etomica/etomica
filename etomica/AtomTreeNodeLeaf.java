package etomica;

/**
 * Leaf node in the tree of atoms.  Differs from group node in having firstChild, lastChild,
 * firstLeafAtom, lastLeafAtom, all given as this node's atom.  Having firstChild
 * and lastChild point to itself is useful in looping through interactions between a leaf
 * and a group.
 */

public final class AtomTreeNodeLeaf extends AtomTreeNode {
    
    /**
     * Linker used to form a list of all leaf atoms in the phase.
     * List is maintained by the speciesMaster node.
     */
    public final AtomLinker leafLinker;
    
    public AtomTreeNodeLeaf(Atom atom, AtomTreeNodeGroup parent) {
        super(atom, parent);
        leafLinker = new AtomLinker(atom);
//        leafCount = (atom.type instanceof AtomType.Wall) ? 0 : 1;
    }
    
    public boolean isLeaf() {return true;}
    
    /**
     * Returns this node's atom.
     */
    public Atom firstLeafAtom() {return atom;}
    
    /**
     * Returns this node's atom.
     */
    public Atom lastLeafAtom() {return atom;}
    
    /**
     * Returns 1.
     */
    public int leafAtomCount() {return 1;}
    
    /**
     * Returns 0.
     */
    public int childAtomCount() {return 0;}

    public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup parent) {
            return new AtomTreeNodeLeaf(atom, parent);
        }
    };
    
}