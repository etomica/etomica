package etomica;

/**
 * Collection of atoms (or other atom groups). 
 *
 * @author David Kofke
 */
public class AtomGroup extends Atom /*implements AtomIteratorBasis */{
    
    protected int childAtomCount;
    protected int leafAtomCount;
    protected boolean resizable;

    /**
     * Constructs an empty atom group with no associated factory.  Normally
     * the new group will be filled with atoms following its construction.
     */
    public AtomGroup(Space space, AtomType.Group type) {
        super(space, type);
        childAtomCount = 0;
        resizable = true;
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
        Atom[] childArray = new Atom[childAtomCount];
        int i=0;
        for(Atom a=firstChildAtom(); a!=null; a=a.nextAtom()) childArray[i++] = a;
        return childArray;
    }
    
    public Atom randomAtom() {return getAtom((int)(Simulation.random.nextDouble()*childAtomCount));}
    
    /**
     * Gets the child atom corresponding to the given index, numbering the first atom as zero
     * and the last atom as Count-1.
     */
    public Atom getAtom(int index) {
        if(index < 0 || index >= childAtomCount) {
            throw new IndexOutOfBoundsException("Index: "+index+
                                                ", Number of childAtoms: "+childAtomCount);
        }
        Atom atom = null;
        if (index < childAtomCount/2) {
            atom = firstChildAtom();
            for (int i = index; i > 0; i--)
                atom = atom.nextAtom();
        } else {
            atom = lastChildAtom();
            for (int i = childAtomCount-1-index; i > 0; i--)
                atom = atom.previousAtom();
        }
        return atom;
    }//end of getAtom
    
    public void setDepth(int d) {
        super.setDepth(d);
        for(Atom atom=firstChildAtom(); atom!=null; atom=atom.nextAtom()) {
            atom.setDepth(d+1);
            if(atom == lastChildAtom()) break;
        }
    }
            
            
    /**
     * Returns the first leaf atom descended from this group.
     */
    public Atom firstLeafAtom() {
        if(childrenAreGroups()) {
            for(Atom a0=firstChildAtom(); a0!=null; a0=a0.nextAtom()) {
                Atom a1 = ((AtomGroup)a0).firstLeafAtom();
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
                Atom a1 = ((AtomGroup)a0).lastLeafAtom();
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
        aNew.setParentGroup(this);
        if(childAtomCount > 0) { //siblings
            aNew.setNextAtom(lastChildAtom().nextAtom());
            lastChildAtom().setNextAtom(aNew);
            aNew.setIndex(lastChildAtom().index()+1);
        }
        else {  //only child
            setFirstAtom(aNew);
            aNew.setNextAtom(null);
            aNew.clearPreviousAtom();
            aNew.setIndex(0);
        }   
        setLastAtom(aNew);
        
        //set up leaf links
        if(aNew instanceof AtomGroup) {
            AtomGroup aNewGroup = (AtomGroup)aNew;//copy to save multiple casts in calls to AtomGroup methods
            Atom lastNewLeaf = aNewGroup.lastLeafAtom();
            if(lastNewLeaf != null) {//also implies firstNewLeaf != null
                lastNewLeaf.setNextAtom(findNextLeafAtom(aNewGroup));
                Atom previous = findPreviousLeafAtom(aNewGroup);
                if(previous != null) previous.setNextAtom(aNewGroup.firstLeafAtom());
            } //no leaf changes anywhere if aNew has no leaf atoms 
        }
        else {//aNew is a leaf, not a group
                //set next in case it wasn't set by aNew.setNextAtom(lastChildAtom.next()) line, above
            if(aNew.nextAtom() == null) aNew.setNextAtom(findNextLeafAtom(aNew.parentGroup()));
               //set previous in case it wasn't set by lastChildAtom.setNextAtom(aNew) line, above
            if(aNew.previousAtom() == null) {
                Atom previous = findPreviousLeafAtom(aNew.parentGroup());
                if(previous != null) previous.setNextAtom(aNew);
            }
        }
        childAtomCount++;
        addAtomNotify(aNew);
    }//end of addAtom


    public void removeAtom(Atom a) {
        if(a == null || !resizable) return;
        //ensure that the given atom is compatible
        if(a.parentGroup() != this) {
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
            Atom first = ((AtomGroup)a).firstLeafAtom();
            Atom last = ((AtomGroup)a).lastLeafAtom();
            next = (last != null) ? last.nextAtom() : null;
            previous = (first != null) ? first.previousAtom() : null;
            if(previous != null) previous.setNextAtom(next);
            else if(next != null) next.clearPreviousAtom();
        }
        
        removeAtomNotify(a);
        a.setParentGroup(null);//must follow notify for SpeciesMaster moleculeCount to be updated
    }//end of removeAtom

    /**
     * Removes all child atoms of this group.  Maintains their
     * internal links, and returns the first child, so they may be recovered 
     * for use elsewhere.
     * Does not remove this group from its parent group.
     */
    public Atom removeAll() {
        if(childAtomCount == 0 || !resizable) return null;
        
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
        childAtomCount = 0;
        return first;
    }//end of removeAll
    
    /**
     * Notifies this atom group that an atom has been added to it or one of its descendants.
     */
    protected void addAtomNotify(Atom atom) {
        leafAtomCount += atom.leafAtomCount();
        if(parentGroup != null) parentGroup.addAtomNotify(atom);
    }
    
    protected void removeAtomNotify(Atom atom) {
        leafAtomCount -= atom.leafAtomCount();
        if(parentGroup != null) parentGroup.removeAtomNotify(atom);
    }
    
    /**
     * Indicates whether the number of child atoms of this group can be changed.
     */
    public boolean isResizable() {return resizable;}
//    public void setResizable(boolean b) {resizable = b;}
    
//    public AtomIterator makeChildAtomIterator() {return new ChildAtomIterator();}
//    public AtomIterator makeLeafAtomIterator() {return new LeafAtomIterator();}
    /**
     * Iterator of the children of this group.
     */
/*    public final AtomIteratorSequential childIterator = new ChildAtomIterator();
    
    public final class ChildAtomIterator extends AtomIteratorSequential {
        public Atom defaultFirstAtom() {return firstChildAtom();}
        public Atom defaultLastAtom() {return lastChildAtom();}
        public boolean contains(Atom a) {return a.parentGroup() == AtomGroup.this;}
    }
    public final class LeafAtomIterator extends AtomIteratorSequential {
        public Atom defaultFirstAtom() {return firstLeafAtom();}
        public Atom defaultLastAtom() {return lastLeafAtom();}
        public boolean contains(Atom a) {return a.isDescendedFrom(AtomGroup.this);}
    }
    */
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public final boolean childrenAreGroups() {
        return firstChildAtom() instanceof AtomGroup;
    }

    /**
    * Chooses a child atom randomly from this group.
    *
    * @return the randomly seleted atom
    */
/*    public Atom randomAtom() {
        int i = (int)(rand.nextDouble()*atomCount);
        Atom a = firstChildAtom;
        for(int j=i; --j>=0; ) {a = a.nextAtom();}
        return a;
    }
*/   

    public final Atom firstChildAtom() {return ((Space.CoordinateGroup)coord).firstAtom();}
    protected final void setFirstAtom(Atom atom) {((Space.CoordinateGroup)coord).setFirstAtom(atom);}
    public final Atom lastChildAtom() {return ((Space.CoordinateGroup)coord).lastAtom();}
    protected final void setLastAtom(Atom atom) {((Space.CoordinateGroup)coord).setLastAtom(atom);}
/*
    //alternative approach that doesn't delegate link structure to coordinates
    private Atom firstChildAtom;
    private Atom lastChildAtom;
    public final Atom firstChildAtom() {return firstChildAtom;}
    protected final void setfirstChildAtom(Atom atom) {firstChildAtom = atom;}
    public final Atom lastChildAtom() {return lastChildAtom;}
    protected final void setlastChildAtom(Atom atom) {lastChildAtom = atom;}
 */   
    /**
     * Searches the atom hierarchy up from (and not including) the given group 
     * to find the closest (leaf) atom in that direction.
     */
    private static Atom findNextLeafAtom(AtomGroup a) {
        AtomGroup parent = a; //search up siblings first
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.nextAtom();
            while(uncle != null) {
                Atom first = uncle.firstLeafAtom();
                if(first != null) return first; //found it
                uncle = (AtomGroup)uncle.nextAtom();
            }
            parent = parent.parentGroup();
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
                Atom last = uncle.lastLeafAtom();
                if(last != null) return last;
                uncle = (AtomGroup)uncle.previousAtom();
            }
            parent = parent.parentGroup();
        }
        return null;
    }//end of findPreviousLeafAtom
    
    private static final IteratorDirective UP = new IteratorDirective(IteratorDirective.UP);
    private static final IteratorDirective DOWN = new IteratorDirective(IteratorDirective.DOWN);
        
}//end of AtomGroup