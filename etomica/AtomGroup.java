package etomica;

/**
 * Collection of atoms (or other atom groups). 
 *
 * @author David Kofke
 */
public class AtomGroup extends Atom /*implements AtomIteratorBasis */{
    
    protected boolean resizable;
//    public int leafAtomCount;
//    public int childAtomCount;

    /**
     * Constructs an empty atom group with no associated factory.  Normally
     * the new group will be filled with atoms following its construction.
     */
    public AtomGroup(Space space, AtomType.Group type) {
        super(space, type);
        resizable = true;
    }
    public AtomGroup(Space space, AtomType.Group type, AtomTreeNode.Factory nodeFactory) {
        super(space, type, nodeFactory);
        resizable = true;
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

/*
    //alternative approach that doesn't delegate link structure to coordinates
    private Atom firstChildAtom;
    private Atom lastChildAtom;
    public final Atom firstChildAtom() {return firstChildAtom;}
    protected final void setfirstChildAtom(Atom atom) {firstChildAtom = atom;}
    public final Atom lastChildAtom() {return lastChildAtom;}
    protected final void setlastChildAtom(Atom atom) {lastChildAtom = atom;}
 */   
    
    private static final IteratorDirective UP = new IteratorDirective(IteratorDirective.UP);
    private static final IteratorDirective DOWN = new IteratorDirective(IteratorDirective.DOWN);
        
}//end of AtomGroup