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
    
    public AtomGroup(AtomGroup parent, int index, AtomFactory factory, 
                        int nChild, Configuration configuration) {
        super(parent, new AtomType.Group(), index);
        this.factory = factory;
        this.configuration = configuration;
        if(nChild > 0) {
            firstAtom = factory.makeAtom(this,0);
            lastAtom = firstAtom;
            for(int i=1; i<nChild; i++) {
                lastAtom.setNextAtom(factory.makeAtom(this,i));
                lastAtom = lastAtom.nextAtom();
                ((Space.CoordinateGroup)coord).addCoordinate(lastAtom.coord);
            }
            factory.configuration().initializeCoordinates(this);
        }
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
    public Atom firstAtom() {
        if(childrenAreGroups()) {
            childIterator.reset(AtomIterator.UP);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).firstAtom();
                if(atom != null) return atom;
            }
            return null;
        }
        else return firstChild;
    }
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastAtom() {
        if(factory.producesAtomGroups()) {
            childIterator.reset(AtomIterator.DOWN);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).lastAtom();
                if(atom != null) return atom;
            }
            return null;
        }
        else return lastChild;
    }

    public void addAtom(Atom a) {
        //ensure that the given atom is compatible
        if(a.parentGroup().atomFactory() != this.factory) {
            System.out.println("Error in AtomGroup.addAtom:  See source code");  //should throw an exception
            return;
        }
        a.setParentGroup(this);
        if(childCount > 0) {//a has siblings
            a.setNextAtom(lastChild.nextAtom());
            lastChild.setNextAtom(a);
            lastChild = a;
        }
        else {  //a is an only child
            firstChild = a;
            lastChild = a;
            a.setNextAtom(null);
            a.clearPreviousAtom();
            //join up the leaf atoms
            Atom atom = findPreviousLastAtom();
            if(previous != null) previous.setNext(a.firstAtom());
            atom = a.lastAtom();
            if(atom != null) atom.setNext(findNextFirstAtom());
        }
        childCount++;
    }


    public void removeAtom(Atom a) {
        //ensure that the given atom is compatible
        if(a.parentGroup() != this) {
            System.out.println("Error in AtomGroup.removeAtom:  See source code");  //should throw an exception
            return;
        }
        a.setParentGroup(this);
        if(childCount > 0) {//a has siblings
            a.setNextAtom(lastChild.nextAtom());
            lastChild.setNextAtom(a);
            lastChild = a;
        }
        else {  //a is an only child
            firstChild = a;
            lastChild = a;
            a.setNextAtom(null);
            a.clearPreviousAtom();
            //join up the leaf atoms
            Atom atom = findPreviousLastAtom();
            if(previous != null) previous.setNext(a.firstAtom());
            atom = a.lastAtom();
            if(atom != null) atom.setNext(findNextFirstAtom());
        }
        childCount++;
    }
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
    * Set the atom following this one, and also sets this one's last 
    * descendent atom to the first descendent atom of the given atom.
    */
    public void setNextAtom(Atom a) {
        //connect to sibling
        super.setNextAtom(a);
        //connect leaf atoms
        Atom lastAtom = lastAtom();
        if(lastAtom != null) {
            if(a instanceof AtomGroup) lastAtom.setNextAtom(a.firstAtom());
            else lastAtom.setNextAtom(a);//handles a being null too
        }
    }
    
    /**
     * Iterator of the children of this group.
     */
    public final AtomIteratorSequential childIterator = new AtomIteratorSequential() {
        public Atom defaultFirstAtom() {return firstChild;}
        public Atom defaultLastAtom() {return lastChild;}
        public boolean contains(Atom a) {return a.parentGroup == AtomGroup.this;}
    };
    
    /**
     * Searches the atom hierarchy up from this group to find the first (leaf) atom.
     * The last leaf atom of this group will have the returned atom as its next atom. 
     */
    private Atom findNextFirstAtom() {
        AtomGroup parent = parentGroup();
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.nextAtom();
            while(uncle != null) {
                Atom first = uncle.firstAtom();
                if(first != null) return first;
                uncle = (AtomGroup)uncle.nextAtom();
            }
            parent = parent.parentGroup();
        }
        return null;
    }//end of findNextFirstAtom
    
    /**
     * Searches the atom hierarchy down from this group to find the last (leaf) atom
     * in that direction.
     * The last first atom of this group will have the returned atom as its previous atom. 
     */
    private Atom findPreviousLastAtom() {
        AtomGroup parent = parentGroup();
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.previousAtom();
            while(uncle != null) {
                Atom last = uncle.lastAtom();
                if(last != null) return last;
                uncle = (AtomGroup)uncle.previousAtom();
            }
            parent = parent.parentGroup();
        }
    }//end of findPreviousLastAtom
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public final boolean childrenAreGroups() {
        return factory.producesAtomGroups();
    }
    
}//end of AtomGroup