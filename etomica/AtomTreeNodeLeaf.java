package etomica;

public final class AtomTreeNodeLeaf implements AtomTreeNode {
    
    public AtomTreeNodeLeaf(Atom atom) {
        setAtom(atom);
    }
    public AtomTreeNodeLeaf() {}
        
    public void setAtom(Atom atom) {
        this.atom = atom;
        leafCount = (atom.type instanceof AtomType.Wall) ? 0 : 1;
    }
    public Atom atom() {return atom;}
        
    public AtomGroup parentGroup() {
        return parentGroup;
   //     return (parentNode != null) ? (AtomGroup)parentNode.atom() : null;
    }
    
    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom parentMolecule() {
        return (parentNode.atom() instanceof SpeciesAgent) ? this.atom : parentNode.parentMolecule();
    }
    
    public void setParentGroup(AtomGroup parent) {
        parentGroup = parent;
        parentNode = parent.node;
        if(parentNode != null) depth = parentNode.depth() + 1;
    }

    public void setDepth(int d) {
        depth = d;
    }
    
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int index() {return atomIndex;}
    public final void setIndex(int i) {atomIndex = i;}

    public final Atom firstChildAtom() {return null;}
    public final Atom lastChildAtom() {return null;}

    public Atom randomAtom() {return null;}
    
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
     /*   AtomGroup ancestor = parentGroup;
        while(ancestor != null) {
            if(ancestor == group) return true;
            ancestor = ancestor.parentGroup();
        }
        return false;
    }*/
        
    public Atom firstLeafAtom() {return atom;}
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastLeafAtom() {return atom;}
    
    public int leafAtomCount() {return leafCount;}
    public int childAtomCount() {return 0;}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        return new Atom[0];
    }
    
            
    /**
     * Returns the first leaf atom descended from this group.
     */
    
    public void addAtom(Atom aNew) {
    }//end of addAtom


    public void removeAtom(Atom a) {
    }//end of removeAtom

    /**
     * Removes all child atoms of this group.  Maintains their
     * internal links, and returns the first child, so they may be recovered 
     * for use elsewhere.
     * Does not remove this group from its parent group.
     */
    public Atom removeAll() {return null;
    }//end of removeAll
    
    /**
     * Notifies this atom group that an atom has been added to it or one of its descendants.
     */
    public void addAtomNotify(Atom atom) {
    }
    
    public void removeAtomNotify(Atom atom) {
    }
    
    private AtomTreeNode parentNode;
    private AtomGroup parentGroup;
    private Atom atom;
    private int leafCount;
    private int depth;
    private int atomIndex;
        
}