package etomica;

/**
 * Leaf node in the tree of atoms.  Differs from group node in having firstChild, lastChild,
 * firstLeafAtom, lastLeafAtom, all given as this node's atom.  Having firstChild
 * and lastChild point to itself is useful in looping through interactions between a leaf
 * and a group.
 */

public final class AtomTreeNodeLeaf implements AtomTreeNode {
    
    /**
     * Linker used to form a list of all leaf atoms in the phase.
     * List is maintained by the speciesMaster node.
     */
    public final AtomLinker leafLinker;
    
    public AtomTreeNodeLeaf(Atom atom) {
        this.atom = atom;
        leafLinker = new AtomLinker(atom);
        leafCount = (atom.type instanceof AtomType.Wall) ? 0 : 1;
    }

    public Atom atom() {return atom;}
        
    public Atom parentGroup() {
        return parentGroup;
   //     return (parentNode != null) ? parentNode.atom() : null;
    }
    
    public AtomTreeNodeGroup parentNode() {
        return parentNode;
    }
    
    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom parentMolecule() {
        return (parentNode.atom() instanceof SpeciesAgent) ? this.atom : parentNode.parentMolecule();
    }
    
    public void setParent(Atom parent) {setParent((AtomTreeNodeGroup)parent.node);}
    
    public void setParent(AtomTreeNodeGroup parent) {
        parentNode = parent;
        parentGroup = (parent != null) ? parent.atom : null;
        if(parentNode != null) {
            depth = parentNode.depth() + 1;
            //parentPhase = parentNode.parentPhase();
        }
        atom.seq.setParentNotify(parent);
    }

    public void setDepth(int d) {
        depth = d;
        if(parentNode != null) parentPhase = parentNode.parentPhase();
    }
    
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int index() {return atomIndex;}
    public final void setIndex(int i) {atomIndex = i;}

    public final Atom firstChildAtom() {return atom;}//null;}
    public final Atom lastChildAtom() {return atom;}//null;}

    public Atom randomAtom() {return atom;}
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public boolean childrenAreGroups() {
        return false;
    }

    /**
     * Gets the child atom corresponding to the given index, numbering the first atom as zero
     * and the last atom as Count-1.
     */
    public Atom getAtom(int index) {
        return null;
    }//end of getAtom
            
    /**
     * Simulation in which this atom resides
     */
    public Simulation parentSimulation() {return parentPhase().parentSimulation();}        
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {
        return parentNode.parentPhase();
    }

    public Species parentSpecies() {return parentNode.parentSpecies();}
    
    public SpeciesAgent parentSpeciesAgent() {return parentNode.parentSpeciesAgent();}

    /**
     * Returns the depth of this atom in the atom hierarchy.  That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public int depth() {return depth;}//return (parentGroup != null) ? parentGroup.depth()+1 : 0;}
    
    public boolean isLeaf() {return true;}
    /**
     * Returns true if the given atom is in the hierarchy of parents of this atom,
     * or if the given atom is this atom.  Returns true, for example, if the given
     * atom is this atom's parent, or its parent's parent, etc.
     */ 
    public boolean isDescendedFrom(Atom group) {
        return (this.atom == group) || (parentNode != null && parentNode.isDescendedFrom(group));
    }
    
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
     * Returns 1.
     */
    public int childAtomCount() {return 1;}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        return new Atom[] {atom};
    }
    
    /**
     * Always throws a RuntimeException because atom cannot be added to leaf.
     */
    public void addAtom(Atom aNew) {
        throw new RuntimeException("Inappropriate call to addAtom in AtomTreeNodeLeaf");
    }//end of addAtom


    /**
     * Always throws a RuntimeException because atom cannot be removed from leaf.
     */
    public void removeAtom(Atom a) {
        throw new RuntimeException("Inappropriate call to addAtom in AtomTreeNodeLeaf");
    }//end of removeAtom

    /**
     * Always throws a RuntimeException.
     */
    public Atom removeAll() {
        throw new RuntimeException("Inappropriate call to addAtomNotify in AtomTreeNodeLeaf");
    }//end of removeAll
    
    /**
     * Always throws a RuntimeException.
     */
    public void addAtomNotify(Atom atom) {
        throw new RuntimeException("Inappropriate call to addAtomNotify in AtomTreeNodeLeaf");
    }
    
    public void removeAtomNotify(Atom atom) {
        throw new RuntimeException("Inappropriate call to removeAtomNotify in AtomTreeNodeLeaf");
    }
    
    private AtomTreeNodeGroup parentNode;
    private Atom parentGroup;
    private Phase parentPhase;
    private final Atom atom;
    private int leafCount;
    private int depth;
    private int atomIndex;
        
    public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
        public AtomTreeNode makeNode(Atom atom) {
            return new AtomTreeNodeLeaf(atom);
        }
    };
    
}