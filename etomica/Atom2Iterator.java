    /**
     * Basic class for iterating over pairs of atoms.
     * Pairs are iterated by collecting pairs yielded by two atom iterators.
     * Different types of pair iterators can be constructed with different choices
     * of the atom iterators.
     */
    public class Atom2Iterator implements java.io.Serializable {
        private AtomPair pair;  //want final, but derived class MP somehow prevents compiler from doing this
        /**
         * The iterators used to generate the sets of atoms
         */
        protected Atom.Iterator ai1, ai2;
        /**
         * A pair action wrapper used to enable the allPairs method
         */
        protected ActionWrapper actionWrapper;   // want final too //wrapper inner class defined below
        protected boolean hasNext;
        /**
         * Flag indicating whether atom1 of pair needs to be updated to point to the same atom that "atom1" in this class points to
         */
        protected boolean needUpdate1; 
        private Atom atom1;
        /**
         * Construct a pair iterator for use in the given phase.  Initial state is hasNext = false
         */
        public Atom2Iterator(Phase p) {
            pair = new AtomPair(p); 
            actionWrapper = new AtomPair.ActionWrapper(pair);
            hasNext = false;
        }
        /**
         * Construct a pair iterator for use in the given phase and which performs the given action.
         * Initial state is hasNext = false
         */
        public Iterator(Phase p, ActionWrapper wrap) {
            pair = new AtomPair(p); 
            actionWrapper = wrap;
            hasNext = false;
        }
        /**
         * Construct a pair iterator for the given phase, using the given atom iterators
         */
        public Iterator(Phase p, Atom.Iterator iter1, Atom.Iterator iter2) {
            pair = new AtomPair(p);
            actionWrapper = new AtomPair.ActionWrapper(pair);
            hasNext = false;
            ai1 = iter1;
            ai2 = iter2;
        }
        public final boolean hasNext() {return hasNext;}
        public void reset() {
            reset(null);
        }
        /**
         * Resets the first atom iterator with the given atom as an argument
         * Then resets the second iterator using the first non-null atom obtained from 
         * the first iterator.
         */
        public void reset(Atom a1) {
            if(a1==null) ai1.reset();
            else ai1.reset(a1);
            do {                  //advance up list of atoms until a pair is found
                if(ai1.hasNext()) {
                    atom1 = ai1.next();
                    ai2.reset(atom1);
                }  
                else {hasNext = false; return;}}   //end of list of atoms
            while(!ai2.hasNext());
            needUpdate1 = true;
            hasNext = true;
        }
        /**
         * Resets the first and second atom iterators using the first and second arguments, respectively.
         */
        public void reset(Atom a1, Atom a2) {
            ai1.reset(a1);
            ai2.reset(a2);
            pair.atom1 = ai1.next();
            needUpdate1 = false;
            hasNext = ai1.hasNext() && ai2.hasNext();
        }
            
        public AtomPair next() {
            if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //ai1 was advanced
            pair.atom2 = ai2.next();
            pair.reset();
            while(!ai2.hasNext()) {     //ai2 is done for this atom1, loop until it is prepared for next
                if(ai1.hasNext()) {     //ai1 has another atom1...
                    atom1 = ai1.next();     //...get it
                    ai2.reset(atom1);       //...reset ai2
                    needUpdate1 = true;     //...flag update of pair.atom1 for next time
                }
                else {hasNext = false; break;} //ai1 has no more; all done with pairs
            }
            return pair;
        }
        
        /**
         * Performs the given action on all pairs returned by this iterator
         */
        public void allPairs(AtomPair.Action act) {  
            reset();
            ai1.reset();  //this shouldn't be here, in general; need to consider it more carefully
            actionWrapper.pairAction = act;
            while(ai1.hasNext()) {
                pair.atom1 = ai1.next();
                ai2.reset(pair.atom1);
                ai2.allAtoms(actionWrapper);
            }
        }
        
        // The following are convenience extensions of AtomPair.Iterator that
        // handle some common iteration needs
        
        /**
         * Iterator for all atom pairs in a phase
         * Default is to do inter and intra pairs; this may be overridden using reset method to do
         * only intermolecular pairs
         * Uses atom iterator and atomPair iterator given by the phase.iteratorFactory class.
         */
         public static final class All extends Iterator {
            public All(Phase p) {
                super(p);
                ai1 = p.iteratorFactory().makeAtomIteratorUp();
                ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
                this.reset();
            }
         }
         
       /**
        * Iterates over pairs formed by given atom and all atoms from other molecules above it in list
        * If given atom is not in phase, it is considered the last atom, and no iterations are performed
        */
         public static final class Up extends Iterator {
            public Up(Phase p) {
                super(p);
                ai1 = new Atom.Iterator.Singlet();
                ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
                this.reset();
            }
            public Up(Phase p, Atom a) {
                super(p);
                ai1 = new Atom.Iterator.Singlet(a);
                ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
                this.reset();
            }
         }
         
       /**
        * Iterates over pairs formed by given atom and all atoms from other molecules below it in list
        * If given atom is not in phase, it is considered the last atom, and iterations are performed
        * over pairs formed from it and all atoms in phase
        */
         public static final class Down extends Iterator {
            public Down(Phase p) {
                super(p);
                ai1 = new Atom.Iterator.Singlet();
                ai2 = p.iteratorFactory().makeAtomIteratorDownNeighbor();
                this.reset();
            }
            public Down(Phase p, Atom a) {
                super(p);
                ai1 = new Atom.Iterator.Singlet(a);
                ai2 = p.iteratorFactory().makeAtomIteratorDownNeighbor();
                this.reset();
            }
         }
        
        /**
         * Iterator for all atoms in a molecule with all atoms in a phase
         * The molecule may or may not be in the phase
         * Intramolecular pairs are not generated
         */
         // Needs to be fixed to handle multi-atom molecules
         public static class MP extends Iterator {
            private boolean upDone;
            private final Atom.Iterator aiUp, aiDown;
            private Atom.Iterator apiCurrent;
            private Molecule molecule;
            public MP(IteratorFactory factory) {
                super(factory.phase());
                aiUp = factory.makeAtomIteratorUpNeighbor();
                aiDown = factory.makeAtomIteratorDownNeighbor();
            }
            public MP(Phase p) {
                super(p);
                aiUp = p.iteratorFactory().makeAtomIteratorUpNeighbor();
                aiDown = p.iteratorFactory().makeAtomIteratorDownNeighbor();
        //        aiUp.setIntra(false);
        //        aiDown.setIntra(false);    (need to implement these methods in Atom.Iterator)
            }
            public MP(Phase p, Molecule m) {
                super(p);
                aiUp = p.iteratorFactory().makeAtomIteratorUpNeighbor();
                aiDown = p.iteratorFactory().makeAtomIteratorDownNeighbor();
        //        aiUp.setIntra(false);
        //        aiDown.setIntra(false);    (need to implement these methods in Atom.Iterator)
                reset(m);
            }
            public void reset() {reset(molecule);}
            public void reset(Molecule m) {
                molecule = m;
                ai1 = m.atomIterator;
                aiUp.reset(m.lastAtom());
                aiDown.reset(m.firstAtom());
                ai2 = aiUp;
                super.reset();
                if(hasNext) {
                    upDone = false;}
                else {
                    ai2 = aiDown;
                    super.reset();
                    upDone = true;
                }
            }
            public AtomPair next() {  //not handling intra/inter in well defined way
                AtomPair dummy = super.next();
                if(!hasNext && !upDone) {
                    ai2 = aiDown;
                    super.reset();
                    upDone = true;
                }
                return pair;   //this was set to proper value when super.next() was called
            }
            public void allPairs(AtomPair.Action act) {
                actionWrapper.pairAction = act;
                while(ai1.hasNext()) {
                    pair.atom1 = ai1.next();
                    aiUp.reset(molecule.lastAtom());
                    aiDown.reset(molecule.firstAtom());
                    aiUp.allAtoms(actionWrapper);
                    aiDown.allAtoms(actionWrapper);
                }
            }
        }  //end of MP 
    }  //end of class Iterator
    
