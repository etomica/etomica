package simulate;

public class AtomPair {
    public Atom atom1, atom2;
    public Space.CoordinatePair cPair;
    public Potential potential;
    public AtomPair() {cPair = null;}
    public AtomPair(Phase phase) {
        cPair = phase.space().makeCoordinatePair(phase.boundary());
    }
    public void reset(Atom a1, Atom a2) {
        atom1 = a1; 
        atom2 = a2;
        reset();
    }
    public void reset() {
        cPair.reset(atom1.coordinate(), atom2.coordinate());
    }
    public void reset(Atom a1, Atom a2, Space.CoordinatePair cp) {
        atom1 = a1;
        atom2 = a2;
        cPair = cp;
    }
    public final double r2() {return cPair.r2();}
    public final double v2() {return cPair.v2();}
    public final double vDotr() {return cPair.vDotr();}
    
    public final Atom atom1() {return atom1;}
    public final Atom atom2() {return atom2;}
    
    /**
     * Interface for classes that generate atom pairs according to various criteria
     */
    public interface Iterator {
        
        public boolean hasNext();
        public AtomPair next();
        public void reset();
        
        public interface M extends Iterator {public void reset(Molecule m);}
        public interface A extends Iterator {  //core iterator, found in each Phase subclass
            public void allDone();
            public void reset(Atom a, boolean intra);
        }
        
        /**
         * Iterator for all atom pairs in a phase
         * Default is to do inter and intra pairs; this may be overridden using reset method to do
         * only intermolecular pairs
         * Uses atom iterator and atomPair iterator given by the phase.iterator class.
         */
         public static class All implements Iterator {
            private final Iterator.A apiUp;
            private final Atom.Iterator atomUp;
            private boolean intra;
            private boolean hasNext;
            private AtomPair thisPair, nextPair;
            public All(Phase p) {
                apiUp = p.iterator.makeAtomPairIteratorUp();
                atomUp = p.iterator.makeAtomIteratorUp();
                reset(true);
            }
            public boolean hasNext() {return hasNext;}
            public AtomPair next() {
                thisPair = nextPair;
                if(apiUp.hasNext()) {nextPair = apiUp.next();}
                else {
                    do {  //advance up list of atoms until one with a pair is found
                        if(atomUp.hasNext()) {apiUp.reset(atomUp.next(),intra);}
                        else {hasNext = false; return thisPair;}}   //end of list of atoms
                    while(!apiUp.hasNext());
                    nextPair = apiUp.next();
                }
                return thisPair;
            }
            public void reset() {reset(intra);}
            public void reset(boolean i) {
                intra = i;
                atomUp.reset();
                do {  //advance up list of atoms until one with a pair is found
                    if(atomUp.hasNext()) {apiUp.reset(atomUp.next(),intra);}
                    else {hasNext = false; return;}}   //end of list of atoms
                while(!apiUp.hasNext());
                nextPair = apiUp.next();
                hasNext = true;
            }
         }
                
        /**
         * Iterator for all atoms in a molecule with all atoms in a phase
         * The molecule may or may not be in the phase
         * Intramolecular pairs are not generated
         */
         // Needs to be fixed to handle multi-atom molecules
        public static class MP implements M {
            private final Iterator.A apiUp, apiDown;
            private Iterator.A apiCurrent;
            private boolean hasNext, upDone;
            private AtomPair nextPair, thisPair;
            public MP(Phase p) {
                apiUp = p.iterator.makeAtomPairIteratorUp();
                apiDown = p.iterator.makeAtomPairIteratorDown();
                hasNext = false;
            }
            public MP(Phase p, Molecule m) {
                apiUp = p.iterator.makeAtomPairIteratorUp();
                apiDown = p.iterator.makeAtomPairIteratorDown();
                reset(m);
            }
            public void reset() {  //needs filling in
                System.out.println("reset in AtomPair.Iterator.MP is not implemented");}  
            public boolean hasNext() {return hasNext;}
            public AtomPair next() {
//                thisPair = nextPair;
/*                if(apiCurrent.hasNext()) {nextPair = apiCurrent.next();}
                else {
                    if(upDone) {hasNext = false;}  //all done
                    else {                         //switch to down iterator
                        apiCurrent = apiDown;
                        upDone = true;
                        hasNext = apiCurrent.hasNext();
                        if(hasNext) nextPair = apiCurrent.next();
                    }
                }*/
                thisPair = apiCurrent.next();
                if(!apiCurrent.hasNext()) {
                    if(upDone) {hasNext = false;}  //all done
                    else {                         //switch to down iterator
                        apiCurrent = apiDown;
                        upDone = true;
                        hasNext = apiCurrent.hasNext();
                    }
                }
                return thisPair;
            }
            public void reset(Molecule m) {
                apiUp.reset(m.lastAtom(),false);     //may be missing atoms on molecule between first and last
                apiDown.reset(m.firstAtom(),false);
                if(apiUp.hasNext()) {
                    apiCurrent = apiUp;
                    upDone = false;}
                else {
                    apiCurrent = apiDown;
                    upDone = true;
                }
                hasNext = apiCurrent.hasNext();
//                if(hasNext) nextPair = apiCurrent.next();
            }
        }
                    
  // Iterates over pairs formed by given atom and all atoms from other molecules above it in list
 // If given atom is not in phase, it is considered the last atom, and no iterations are performed
        public static final class Up implements Iterator.A {
            private final AtomPair pair;
            private final Phase phase;
            private boolean hasNext;
            private Atom nextAtom;
            public Up(Phase p) {phase = p; pair = new AtomPair(p); hasNext = false;}
            public Up(Phase p, Atom a) {phase = p; pair = new AtomPair(p); reset(a,true);}
            public void reset() {System.out.println("error in APIup");}
            public void allDone() {hasNext = false;}
            public boolean hasNext() {return hasNext;}
            public void reset(Atom a, boolean intra) {
                if(a == null || a.parentPhase() != phase) {hasNext = false; return;}
                pair.atom1 = a;
                nextAtom = intra ? a.nextAtom() : a.nextMoleculeFirstAtom();
                hasNext = (nextAtom != null);
            }
            public AtomPair next() {
                pair.atom2 = nextAtom;
                pair.reset();
                nextAtom = nextAtom.nextAtom();
                hasNext = (nextAtom != null);
                return pair;
            }
        } //end of AtomPair.Iterator.Up
 
    /**
     * Iterates over pairs formed by given atom and all atoms from other molecules below it in list
     * If given atom is not in phase, it is considered the last atom, and iterations are performed over
     * all atoms in phase
     */
    
        public static final class Down implements Iterator.A {
            private final AtomPair pair;
            private final Phase phase;
            private boolean hasNext;
            private Atom nextAtom;
            public Down(Phase p) {phase = p; pair = new AtomPair(p); hasNext = false;}
            public Down(Phase p, Atom a) {phase = p; pair = new AtomPair(p); reset(a,true);}
            public boolean hasNext() {return hasNext;}
            public void reset() {System.out.println("error in APIdown");}
            public void allDone() {hasNext = false;}
            public void reset(Atom a, boolean intra) {
                if(a == null) {hasNext = false; return;}
                pair.atom1 = a;
                if(a.parentPhase() == phase) {
                    nextAtom = intra ? a.previousAtom() : a.previousMoleculeLastAtom();}
                else {
                    nextAtom = phase.lastAtom();
                }
                hasNext = (nextAtom != null);
            }
            public AtomPair next() {
                pair.atom2 = nextAtom;
                pair.reset();
                nextAtom = nextAtom.previousAtom();
                hasNext = (nextAtom != null);
                return pair;
            }
        } //end of AtomPair.Iterator.Down
        
    }  //end of interface Iterator
}  //end of  AtomPair