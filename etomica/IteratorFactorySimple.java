package etomica;

public class IteratorFactorySimple implements IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
    
    public AtomIterator makeGroupIteratorSimple() {return new AtomIteratorListSimple();}
    
    public AtomIterator makeAtomIterator() {return new Iterator();}
//    public AtomIterator makeAtomIterator() {return new AtomIteratorSequential();}

    public AtomSequencer makeAtomSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
        
    public AtomIterator makeIntragroupIterator() {return new Iterator();}
    public AtomIterator makeIntergroupIterator() {return new Iterator();}
 //   public AtomIterator makeNeighborIterator() {return new AtomIteratorSequential();}
    
    public AtomSequencer makeNeighborSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
    
    public Class atomSequencerClass() {return AtomSequencerSimple.class;}
    
    public Class neighborSequencerClass() {return AtomSequencerSimple.class;}
    
    public AtomSequencer.Factory atomSequencerFactory() {return AtomSequencerSimple.FACTORY;}
    
    public AtomSequencer.Factory neighborSequencerFactory() {return AtomSequencerSimple.FACTORY;}
   
/////////////////////////////////////////////////////////////////////////////////////////////
    /**
    * Iterates over all the children of a given atom group.
    *
    * @author David Kofke
    */

    //won't work correctly if tabs are in sequence
    public static final class Iterator implements AtomIterator {
        
        /**
        * Indicates if another iterate is forthcoming.
        */
        public boolean hasNext() {return next != null;}
        
        /**
        * True if the parent group of the given atom is the current basis for the iterator.
        * False otherwise, or if atom or basis is null.
        */
        public boolean contains(Atom atom) {
            return atom != null && atom.node.isDescendedFrom(basis);
        }
        
        /**
        * Does reset if relation (preceeds) of atom in iterator directive and basis of this
        * iterator are compatible with the direction in iterator directive.  Otherwise
        * resets to hasNext false.
        */
        public Atom reset(IteratorDirective id) {
            switch(id.atomCount()) {
                case 0: return reset();
                case 1: upListNow = id.direction().doUp();
                        doGoDown = id.direction().doDown();
                        return reset(id.atom1());
                default: throw new RuntimeException("Illegal atomCount in iteratorDirective in IteratorFactorySimple.reset(IteratorDirective)");
            }
        }
        
        private Atom reset(Atom atom) {
            if(atom == null) {
                next = null;
                return null;
            }
            
            if(atom == basis) return reset();
            
            if(!this.contains(atom)) {
                boolean before = atom.seq.preceeds(basis);
                if(before && upListNow) return reset();
                else if(!before && doGoDown) {
                    upListNow = false;
                    doGoDown = false;
                    next = basis.node.lastChildAtom();
                    last = basis.node.firstChildAtom();
                    return next;
                }
                else {
                    next = null;
                    return null;
                }
            }
            
            ref = atom;
            next = atom;
            last = upListNow ? basis.node.lastChildAtom() : basis.node.firstChildAtom();
            next();
            return next;
        }
        
        /**
        * Ignored.
        */
        public void setAsNeighbor(boolean b) {}
        
        /**
        * Resets iterator to loop from first child to last child of basis.
        */
        public Atom reset() {
            if(basis == null) {next = null; return null;}
            upListNow = true;
            doGoDown = false;
            next = basis.node.firstChildAtom();
            last = basis.node.lastChildAtom();
            return next;
        }
        
        /**
        * Returns the next atom in the iteration.
        */
        public Atom next() {
            Atom atom = next;
            if(next != last) next = upListNow ? next.seq.nextAtom() : next.seq.previousAtom();
            else if(doGoDown) {
                upListNow = false;
                doGoDown = false;
                next = ref.seq.previousAtom();
                last = basis.node.firstChildAtom();
            }
            else next = null;
            return atom;
        }        
        
        /**
        * Performs given action for each child atom of basis.
        */
        public void allAtoms(AtomAction act) {;
            if(next == null) return;
            if(upListNow) {
                for(Atom atom = next; atom != null; atom=atom.seq.nextAtom()) {
                    act.actionPerformed(atom);
                    if(atom == last) break;
                }
                if(ref == null) {next = null; return;}
                next = ref.seq.previousAtom();
                last = basis.node.firstChildAtom();
            }
            if(doGoDown && next != null) {
                for(Atom atom = next; atom != null; atom=atom.seq.previousAtom()) {
                    act.actionPerformed(atom);
                    if(atom == last) break;
                }
            }
            next = null;
        }
            
        /**
        * Sets the given atom as the basis, so that child atoms of the
        * given atom will be returned upon iteration.  If given atom is
        * a leaf atom, no iterates are given.
        */
        public void setBasis(Atom atom) {
            basis = /*(atom == null || atom.node.isLeaf()) ? null :*/ atom;
        }
        
        /**
        * Returns the current iteration basis.
        */
        public Atom getBasis() {return basis;}
        
        /**
        * The number of atoms returned on a full iteration, using the current basis.
        */
        public int size() {return (basis != null) ? basis.node.childAtomCount() : 0;}   

        private Atom basis;
        private Atom next;
        private Atom last;
        private Atom ref;
        private boolean upListNow, doGoDown;

    }//end of Iterator
/////////////////////////////////////////////////////////////////////////////////////////////

}