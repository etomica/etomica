package etomica;

/**
 * Collection of atoms (or other atom groups). 
 *
 * @author David Kofke
 */
public class AtomGroup extends Atom {
    
    protected final AtomFactory factory;
    private Configuration initializer;
    protected int childCount;
    protected final boolean resizable;
    private static final IteratorDirective UP = new IteratorDirective(IteratorDirective.UP);
    private static final IteratorDirective DOWN = new IteratorDirective(IteratorDirective.DOWN);
    
    //used to pause when adding/removing atoms
    private transient Integrator integrator;
    
    protected AtomGroup(AtomGroup parent, int index, AtomType type) {
        super(parent, index, type);
        childCount = 0;
        resizable = true;
        factory = null;
    }
    public AtomGroup(AtomGroup parent, int index, AtomType type, AtomFactory factory, 
                        int nChild, Configuration initializer, boolean resizable) {
        super(parent, index, type);
        this.factory = factory;
        this.initializer = initializer;
        for(int i=0; i<nChild; i++) {
            addAtom(factory.makeAtom(this,i));
        }
        if(initializer != null) initializer.initializeCoordinates(this); //be sure this handles case of no atoms
        childCount = nChild;
        this.resizable = resizable;
    }
    
    public AtomFactory atomFactory() {return factory;}
    
    public void setConfiguration(Configuration c) {}
    
    //to be completed
    public int atomCount() {return -1;}
    public int childCount() {return childCount;}
    
    /**
     * Returns the first leaf atom descended from this group.
     */
    public Atom firstLeafAtom() {
        if(childrenAreGroups()) {
            childIterator.reset(UP);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).firstLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        }
        else return firstChild();
    }
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastLeafAtom() {
        if(childrenAreGroups()) {
            childIterator.reset(DOWN);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).lastLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        }
        else return lastChild();
    }
    
    public Atom addNewAtom() {
        Atom aNew = atomFactory().makeAtom(this, childCount+1);
        addAtom(aNew);
        return aNew;
    }

    public void addAtom(Atom aNew) {
        if(!resizable) return;
        //ensure that the given atom is compatible
        if(atomFactory().vetoAddition(aNew)) { //also checks for a null
            System.out.println("Error in AtomGroup.addAtom:  See source code");  //should throw an exception
            return;
        }
        //set up links within this group
        aNew.setParentGroup(this);
        if(childCount > 0) { //siblings
            aNew.setNextAtom(lastChild().nextAtom());
            lastChild().setNextAtom(aNew);
            aNew.setIndex(lastChild().index()+1);
        }
        else {  //only child
            setFirstChild(aNew);
            aNew.setNextAtom(null);
            aNew.clearPreviousAtom();
            aNew.setIndex(0);
        }   
        setLastChild(aNew);
        
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
                //set next in case it wasn't set by aNew.setNextAtom(lastChild.next()) line, above
            if(aNew.nextAtom() == null) aNew.setNextAtom(findNextLeafAtom(aNew.parentGroup()));
               //set previous in case it wasn't set by lastChild.setNextAtom(aNew) line, above
            if(aNew.previousAtom() == null) {
                Atom previous = findPreviousLeafAtom(aNew.parentGroup());
                if(previous != null) previous.setNextAtom(aNew);
            }
        }
        childCount++;
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
        if(a == firstChild()) {
            if(childCount == 1) setFirstChild(null);
            else setFirstChild(next);
        }
        if(a == lastChild()) {
            if(childCount == 1) setLastChild(null);
            else setLastChild(previous);
        }
        if(previous != null) previous.setNextAtom(next);
        else if(next != null) next.clearPreviousAtom();
        childCount--;
        
        a.setNextAtom(null);
        a.clearPreviousAtom();        
        a.setParentGroup(null);
        
        //update leaf-atom links
        if(a instanceof AtomGroup) {
            Atom first = ((AtomGroup)a).firstLeafAtom();
            Atom last = ((AtomGroup)a).lastLeafAtom();
            next = (last != null) ? last.nextAtom() : null;
            previous = (first != null) ? first.previousAtom() : null;
            if(previous != null) previous.setNextAtom(next);
            else if(next != null) next.clearPreviousAtom();
        }
        
        factory.reservoir().addAtom(a);
    }//end of removeAtom

    /**
     * Removes all child atoms of this group.  
     * Does not remove this group from its parent group.
     */
    public void removeAll() {
        if(childCount == 0 || !resizable) return;
        
        Atom first = firstLeafAtom();
        Atom last = lastLeafAtom();
        Atom next = (last != null) ? last.nextAtom() : null;
        Atom previous = (first != null) ? first.previousAtom() : null;
        if(previous != null) previous.setNextAtom(next);
        else if(next != null) next.clearPreviousAtom();
        
        setFirstChild(null);
        setLastChild(null);
        childCount = 0;
        
    }//end of removeAll
    
    /**
     * Indicates whether the number of child atoms of this group can be changed.
     */
    public boolean isResizable() {return resizable;}
//    public void setResizable(boolean b) {resizable = b;}
    
    /**
     * Iterator of the children of this group.
     */
    public final AtomIteratorSequential childIterator = new AtomIteratorSequential() {
        public Atom defaultFirstAtom() {return firstChild();}
        public Atom defaultLastAtom() {return lastChild();}
        public boolean contains(Atom a) {return a.parentGroup() == AtomGroup.this;}
    };
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public final boolean childrenAreGroups() {
        return factory.producesAtomGroups();
    }

    
    
    /**
    * Sets the number of molecules for this species.  Makes the given number
    * of new molecules, linked-list orders and initializes them.
    * Any previously existing molecules for this species in this phase are abandoned
    * Any links to molecules of next or previous species are maintained.
    * Takes no action at all if the new number of molecules equals the existing number
    *
    * @param n  the new number of molecules for this species
    * @see #makeMolecule
    * @see #deleteMolecule
    * @see #addMolecule
    */
    public void setNAtoms(int n) {
        if(!resizable) return;
        boolean wasPaused = pauseIntegrator();
        
        if(n <= 0) removeAll();
        else if(n > childCount) {
            for(int i=childCount; i<n; i++) addNewAtom();
        }
        else if(n < childCount) {
            for(int i=childCount; i>n; i--) removeAtom(lastChild());
        }
        
        //reconsider this
        parentPhase().configuration.initializeCoordinates(this);
        parentPhase().iteratorFactory().reset();
        
        unpauseIntegrator(wasPaused);
    }
        
    /**
     * Same as setNAtoms, but takes a boolean argument that can indicate that all
     * new molecules should be made.
     */
    public void setNAtoms(int n, boolean forceRebuild) {
        if(!resizable) return;
        boolean wasPaused = pauseIntegrator();
        if(forceRebuild) removeAll();
        setNAtoms(n);
        unpauseIntegrator(wasPaused);
    }
    
    private boolean pauseIntegrator() {
        Phase phase = parentPhase();
        integrator = (phase != null) ? phase.integrator() : null;
        boolean wasPaused = true;
        if(integrator != null) {
            wasPaused = integrator.isPaused();//record pause state of integrator
            if(!wasPaused) {
                integrator.pause();
                while(!integrator.isPaused()) {}
            }
        }
        return wasPaused;
    }
    
    private void unpauseIntegrator(boolean wasPaused) {
        if(integrator != null) {
            if(integrator.isInitialized()) integrator.initialize();//reinitialize only if initialized already
            if(!wasPaused) integrator.unPause();//resume if was not paused originally
        }
    }
              
    /**
    * Chooses a child atom randomly from this group.
    *
    * @return the randomly seleted atom
    */
/*    public Atom randomAtom() {
        int i = (int)(rand.nextDouble()*childCount);
        Atom a = firstChild;
        for(int j=i; --j>=0; ) {a = a.nextAtom();}
        return a;
    }
*/   

    public final Atom firstChild() {return ((Space.CoordinateGroup)coord).firstChild();}
    protected final void setFirstChild(Atom atom) {((Space.CoordinateGroup)coord).setFirstChild(atom);}
    public final Atom lastChild() {return ((Space.CoordinateGroup)coord).lastChild();}
    protected final void setLastChild(Atom atom) {((Space.CoordinateGroup)coord).setLastChild(atom);}
/*
    //alternative approach that doesn't delegate link structure to coordinates
    private Atom firstChild;
    private Atom lastChild;
    public final Atom firstChild() {return firstChild;}
    protected final void setFirstChild(Atom atom) {firstChild = atom;}
    public final Atom lastChild() {return lastChild;}
    protected final void setLastChild(Atom atom) {lastChild = atom;}
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
    
}//end of AtomGroup