package etomica;

public class IteratorFactorySimple implements IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
    
 //   public AtomIterator makeGroupIteratorSimple() {return new AtomIteratorListSimple();}
    
    public AtomIterator makeGroupIteratorSequential() {return new Iterator();}
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
            return atom != null && atom.node.isDescendedFrom(basisNode.atom);
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
            
            if(atom == basisNode.atom) return reset();
            
            if(!this.contains(atom)) {
                boolean before = atom.seq.preceeds(basisNode.atom);
                if(before && upListNow) return reset();
                else if(!before && doGoDown) {
                    upListNow = false;
                    doGoDown = false;
                    next = basisNode.lastChildAtom();
                    last = basisNode.firstChildAtom();
                    return next;
                }
                else {
                    next = null;
                    return null;
                }
            }
            
            ref = atom;
            next = atom;
            last = upListNow ? basisNode.lastChildAtom() : basisNode.firstChildAtom();
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
            if(basisNode == null) {next = null; return null;}
            upListNow = true;
            doGoDown = false;
            next = basisNode.firstChildAtom();
            last = basisNode.lastChildAtom();
            return next;
        }
        
        /**
        * Returns the next atom in the iteration.
        */
        public Atom next() {
            Atom atom = next;
            if(next != last) next = upListNow ? next.seq.next.atom : next.seq.previous.atom;
            else if(doGoDown) {
                upListNow = false;
                doGoDown = false;
                next = ref.seq.previous.atom;
                last = basisNode.firstChildAtom();
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
                for(Atom atom = next; atom != null; atom=atom.seq.next.atom) {
                    act.actionPerformed(atom);
                    if(atom == last) break;
                }
                if(ref == null) {next = null; return;}
                next = ref.seq.previous.atom;
                last = basisNode.firstChildAtom();
            }
            if(doGoDown && next != null) {
                for(Atom atom = next; atom != null; atom=atom.seq.previous.atom) {
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
            basisNode = (atom == null || atom.node.isLeaf()) ? null : (AtomTreeNodeGroup)atom.node;
        }
        
        /**
        * Returns the current iteration basis.
        */
        public Atom getBasis() {return basisNode.atom;}
        
        /**
        * The number of atoms returned on a full iteration, using the current basis.
        */
        public int size() {return (basisNode != null) ? basisNode.childAtomCount() : 0;}   

        private AtomTreeNodeGroup basisNode;
        private Atom next;
        private Atom last;
        private Atom ref;
        private boolean upListNow, doGoDown;

    }//end of Iterator
/////////////////////////////////////////////////////////////////////////////////////////////

}