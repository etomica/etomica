package etomica;

/**
 * Implementation of AtomTreeNode that relies on recursion to
 * evaluate most elements of hierarchy relative to the atom.
 */
 
public class AtomTreeNodeRecursive implements AtomTreeNode {
   
    public AtomTreeNodeRecursive(Atom atom) {
        setAtom(atom);
    }
    public AtomTreeNodeRecursive() {}
        
    public void setAtom(Atom atom) {
        this.atom = atom;
    }
    public Atom atom() {return atom;}
        
    public AtomGroup parentGroup() {
        return parentGroup;
 //       return (parentNode != null) ? (AtomGroup)parentNode.atom() : null;
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
        for(Atom atom=firstChildAtom(); atom!=null; atom=atom.nextAtom()) {
            atom.node.setDepth(d+1);
            if(atom == lastChildAtom()) break;
        }
    }
    
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int index() {return atomIndex;}
    public final void setIndex(int i) {atomIndex = i;}
 
    public final Atom firstChildAtom() {return ((Space.CoordinateGroup)atom.coord).firstAtom();}
    public final Atom lastChildAtom() {return ((Space.CoordinateGroup)atom.coord).lastAtom();}

    public Atom randomAtom() {return getAtom((int)(Simulation.random.nextDouble()*childAtomCount()));}
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public boolean childrenAreGroups() {
        return !firstChildAtom().node.isLeaf();
    }

    /**
     * Gets the child atom corresponding to the given index, numbering the first atom as zero
     * and the last atom as Count-1.
     */
    public Atom getAtom(int index) {
        if(index < 0 || index >= childAtomCount()) {
            throw new IndexOutOfBoundsException("Index: "+index+
                                                ", Number of childAtoms: "+childAtomCount());
        }
        Atom atom = null;
        if (index < childAtomCount()/2) {
            atom = firstChildAtom();
            for (int i = index; i > 0; i--)
                atom = atom.nextAtom();
        } else {
            atom = lastChildAtom();
            for (int i = childAtomCount()-1-index; i > 0; i--)
                atom = atom.previousAtom();
        }
        return atom;
    }//end of getAtom
            
    /**
     * Simulation in which this atom resides
     */
    public Simulation parentSimulation() {return parentPhase().parentSimulation();}        
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {return parentNode.parentPhase();}

    public Species parentSpecies() {return parentNode.parentSpecies();}
    
    public SpeciesAgent parentSpeciesAgent() {return parentNode.parentSpeciesAgent();}

    /**
     * Returns the depth of this atom in the atom hierarchy.  That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public int depth() {return depth;}//return (parentGroup != null) ? parentGroup.depth()+1 : 0;}
    
    public boolean isLeaf() {return false;}
    
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
        
    public Atom firstLeafAtom() {
        if(childrenAreGroups()) {
            for(Atom a0=firstChildAtom(); a0!=null; a0=a0.nextAtom()) {
                Atom a1 = a0.node.firstLeafAtom();
                if(a1 != null) return a1;
            }
            return null;
        }
          /* using iterator for loop seems to cause problem with return of null in AtomPairIterator inner loop even though it reports hasNext=true  
            childIterator.reset(UP);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).firstLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        }*/
        else return firstChildAtom();
    }
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastLeafAtom() {
        if(childrenAreGroups()) {
            for(Atom a0=lastChildAtom(); a0!=null; a0=a0.previousAtom()) {
                Atom a1 = ((AtomGroup)a0).node.lastLeafAtom();
                if(a1 != null) return a1;
            }
            return null;
        }
          /* using iterator for loop seems to cause problem with return of null in AtomPairIterator inner loop even though it reports hasNext=true  
            childIterator.reset(DOWN);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).lastLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        } */
        else return lastChildAtom();
    }
    
    public int leafAtomCount() {return leafAtomCount;}
    public int childAtomCount() {return childAtomCount;}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        Atom[] childArray = new Atom[childAtomCount()];
        int i=0;
        for(Atom a=firstChildAtom(); a!=null; a=a.nextAtom()) childArray[i++] = a;
        return childArray;
    }
    
            
    /**
     * Returns the first leaf atom descended from this group.
     */
    
    public void addAtom(Atom aNew) {
        if(!resizable) return; //should define an exception
        //ensure that the given atom is compatible
        //creator() is not defined for speciesAgent, so this creates a problem when adding
        //new molecule; commented out vetoAddition call until this is fixed or abandoned.
        if(aNew == null /* || creator().vetoAddition(aNew)*/) { //should ensure group/atom consistency with other children
            System.out.println("Error in AtomGroup.addAtom:  See source code");  //should throw an exception
            return;
        }
        //set up links within this group
        aNew.node.setParentGroup((AtomGroup)atom);
        if(childAtomCount > 0) { //siblings
            aNew.setNextAtom(lastChildAtom().nextAtom());
            lastChildAtom().setNextAtom(aNew);
            aNew.node.setIndex(lastChildAtom().node.index()+1);
        }
        else {  //only child
            setFirstAtom(aNew);
            aNew.setNextAtom(null);
            aNew.clearPreviousAtom();
            aNew.node.setIndex(0);
        }   
        setLastAtom(aNew);
        
        //set up leaf links
        if(aNew instanceof AtomGroup) {
            AtomGroup aNewGroup = (AtomGroup)aNew;//copy to save multiple casts in calls to AtomGroup methods
            Atom lastNewLeaf = aNewGroup.node.lastLeafAtom();
            if(lastNewLeaf != null) {//also implies firstNewLeaf != null
                lastNewLeaf.setNextAtom(findNextLeafAtom(aNewGroup));
                Atom previous = findPreviousLeafAtom(aNewGroup);
                if(previous != null) previous.setNextAtom(aNewGroup.node.firstLeafAtom());
            } //no leaf changes anywhere if aNew has no leaf atoms 
        }
        else {//aNew is a leaf, not a group
                //set next in case it wasn't set by aNew.setNextAtom(lastChildAtom.next()) line, above
            if(aNew.nextAtom() == null) aNew.setNextAtom(findNextLeafAtom(aNew.node.parentGroup()));
               //set previous in case it wasn't set by lastChildAtom.setNextAtom(aNew) line, above
            if(aNew.previousAtom() == null) {
                Atom previous = findPreviousLeafAtom(aNew.node.parentGroup());
                if(previous != null) previous.setNextAtom(aNew);
            }
        }
        childAtomCount++;
        addAtomNotify(aNew);
    }//end of addAtom


    public void removeAtom(Atom a) {
        if(a == null || !resizable) return;
        //ensure that the given atom is compatible
        if(a.node.parentGroup() != atom) {
            System.out.println("Error in AtomGroup.removeAtom:  See source code");  //should throw an exception
            return;
        }
        Atom next = a.nextAtom();
        Atom previous = a.previousAtom();
        
        //update links within this group
        if(a == firstChildAtom()) {
            if(childAtomCount == 1) setFirstAtom(null);
            else setFirstAtom(next);
        }
        if(a == lastChildAtom()) {
            if(childAtomCount == 1) setLastAtom(null);
            else setLastAtom(previous);
        }
        if(previous != null) previous.setNextAtom(next);
        else if(next != null) next.clearPreviousAtom();
        childAtomCount--;
        
        a.setNextAtom(null);
        a.clearPreviousAtom();        
        
        //update leaf-atom links
        if(a instanceof AtomGroup) {
            Atom first = ((AtomGroup)a).node.firstLeafAtom();
            Atom last = ((AtomGroup)a).node.lastLeafAtom();
            next = (last != null) ? last.nextAtom() : null;
            previous = (first != null) ? first.previousAtom() : null;
            if(previous != null) previous.setNextAtom(next);
            else if(next != null) next.clearPreviousAtom();
        }
        
        removeAtomNotify(a);
        a.node.setParentGroup(null);//must follow notify for SpeciesMaster moleculeCount to be updated
    }//end of removeAtom

    /**
     * Removes all child atoms of this group.  Maintains their
     * internal links, and returns the first child, so they may be recovered 
     * for use elsewhere.
     * Does not remove this group from its parent group.
     */
    public Atom removeAll() {
        if(childAtomCount() == 0 || !resizable) return null;
        
        //disconnect leaf atoms from their their links to adjacent groups
        Atom first = firstLeafAtom();
        Atom last = lastLeafAtom();
        Atom next = (last != null) ? last.nextAtom() : null;
        Atom previous = (first != null) ? first.previousAtom() : null;
        if(previous != null) previous.setNextAtom(next);
        else if(next != null) next.clearPreviousAtom();
        
        //save first child for return, then erase links to them
        first = firstChildAtom();
        setFirstAtom(null);
        setLastAtom(null);
        return first;
    }//end of removeAll
    
    /**
     * Notifies this atom group that an atom has been added to it or one of its descendants.
     */
    public void addAtomNotify(Atom atom) {
        leafAtomCount += atom.node.leafAtomCount();
        if(parentGroup() != null) parentGroup().node.addAtomNotify(atom);
    }
    
    public void removeAtomNotify(Atom atom) {
        leafAtomCount -= atom.node.leafAtomCount();
        if(parentGroup() != null) parentGroup().node.removeAtomNotify(atom);
    }
    

    private final void setFirstAtom(Atom a) {((Space.CoordinateGroup)atom.coord).setFirstAtom(a);}
    private final void setLastAtom(Atom a) {((Space.CoordinateGroup)atom.coord).setLastAtom(a);}
 
    /**
     * Searches the atom hierarchy up from (and not including) the given group 
     * to find the closest (leaf) atom in that direction.
     */
    private static Atom findNextLeafAtom(AtomGroup a) {
        AtomGroup parent = a; //search up siblings first
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.nextAtom();
            while(uncle != null) {
                Atom first = uncle.node.firstLeafAtom();
                if(first != null) return first; //found it
                uncle = (AtomGroup)uncle.nextAtom();
            }
            parent = parent.node.parentGroup();
        }
        return null;
    }//end of findNextLeafAtom
    
    /**
     * Searches the atom hierarchy down from (and not including) the given group 
     * to find the closest (leaf) atom in that direction.
     * The first leaf atom of the given group will have the returned atom as its previous atom. 
     */
    private static Atom findPreviousLeafAtom(AtomGroup a) {
        AtomGroup parent = a;
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.previousAtom();
            while(uncle != null) {
                Atom last = uncle.node.lastLeafAtom();
                if(last != null) return last;
                uncle = (AtomGroup)uncle.previousAtom();
            }
            parent = parent.node.parentGroup();
        }
        return null;
    }//end of findPreviousLeafAtom

    protected Atom atom;
    protected int leafAtomCount;
    
    protected int childAtomCount;
    protected int depth;
    protected int atomIndex;
    private AtomTreeNode parentNode;
    private AtomGroup parentGroup;
    
    private boolean resizable = true;
    
    
}