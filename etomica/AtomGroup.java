package etomica;

/**
 * Collection of atoms (or other atom groups). 
 *
 * @author David Kofke
 */
public class AtomGroup extends Atom implements java.io.Serializable {
    
    private Atom firstChild;
    private Atom lastChild;
    private final AtomFactory factory;
    private int childCount;
    private boolean resizable = true;
    private static final IteratorDirective UP = new IteratorDirective(IteratorDirective.UP);
    private static final IteratorDirective DOWN = new IteratorDirective(IteratorDirective.DOWN);
    
    public AtomGroup(AtomGroup parent, int index, AtomFactory factory, 
                        int nChild, Configuration configuration) {
        super(parent, new AtomType.Group(), index);
        this.factory = factory;
        this.configuration = configuration;
        for(int i=0; i<nChild; i++) {
            addAtom(factory.makeAtom(this,i));
            ((Space.CoordinateGroup)coord).addCoordinate(lastAtom.coord);
        }
        factory.configuration().initializeCoordinates(this); //be sure this handles case of no atoms
        childCount = nChild;
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
        else return firstChild;
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
        else return lastChild;
    }

    public void addAtom(Atom aNew) {
        //ensure that the given atom is compatible
        if(atomFactory().vetoAddition(aNew)) { //also checks for a null
            System.out.println("Error in AtomGroup.addAtom:  See source code");  //should throw an exception
            return;
        }
        //set up links within this group
        aNew.setParentGroup(this);
        if(childCount > 0) { //siblings
            aNew.setNextAtom(lastChild.nextAtom());
            lastChild.setNextAtom(aNew);
            aNew.setIndex(lastChild.index()+1);
        }
        else {  //only child
            firstChild = aNew;
            aNew.setNextAtom(null);
            aNew.clearPreviousAtom();
            aNew.setIndex(0);
        }   
        lastChild = aNew;
        
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
        if(a == null) return;
        //ensure that the given atom is compatible
        if(a.parentGroup() != this) {
            System.out.println("Error in AtomGroup.removeAtom:  See source code");  //should throw an exception
            return;
        }
        Atom next = a.nextAtom();
        Atom previous = a.previousAtom();
        
        //update links within this group
        if(a == firstChild) {
            if(childCount == 1) firstChild = null;
            else firstChild = next;
        }
        if(a == lastChild) {
            if(childCount == 1) lastChild = null;
            else lastChild = previous;
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
    }//end of removeAtom
/*      
    public boolean isResizable() {return resizable;}
    public void setResizable(boolean b) {
        resizable = b;
        if(resizable) parentGroup.setResizable(true);//if this is resizable, so must its parents
        else if(factory.producesAtomGroups()) {//if this isn't resizable, neither can its children
            childIterator.reset();
            while(childIterator.hasNext()) {
                ((AtomGroup)childIterator.next()).setResizable(false);
            }
        }
    }//end setResizable
   */
   
    /**
     * Iterator of the children of this group.
     */
    public final AtomIteratorSequential childIterator = new AtomIteratorSequential() {
        public Atom defaultFirstAtom() {return firstChild;}
        public Atom defaultLastAtom() {return lastChild;}
        public boolean contains(Atom a) {return a.parentGroup() == AtomGroup.this;}
    };
    
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
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public final boolean childrenAreGroups() {
        return factory.producesAtomGroups();
    }
    
}//end of AtomGroup