package simulate;

public interface AtomPairIterator {
    
    public boolean hasNext();
    public AtomPair next();
    
    public interface SS extends AtomPairIterator {public void reset(Species s1, Species s2);}
    public interface S extends AtomPairIterator {public void reset(Species s);}
    public interface MM extends AtomPairIterator {public void reset(Molecule m1, Molecule m2);}
    public interface M extends AtomPairIterator {public void reset(Molecule m);}
    public interface A extends AtomPairIterator {
        public void reset(Atom a1, Atom a2, Atom a3);
        public void reset(Atom a1, Atom a2, Atom a3, Atom a4);
        public void allDone();
    }

    // Iterator for all (intramolecular) atom pairs in a single molecule
    public static class AllMoleculeIntraPairs implements AtomPairIterator.M {
        private final AtomPairIterator.A iterator;
        public AllMoleculeIntraPairs(PhaseSpace ps) {iterator = ps.makePairIteratorHalf();} //constructor
        public AllMoleculeIntraPairs(PhaseSpace ps, Molecule m) {  //constructor
            iterator = ps.makePairIteratorHalf(); 
            reset(m);
        } 
        public boolean hasNext() {return iterator.hasNext();}
        public AtomPair next() {return iterator.next();}
        public final void reset(Molecule m) {
            if(m == null || m.nAtoms < 2) {iterator.allDone(); return;}
            Atom oF = m.firstAtom();
            Atom iF = oF.nextAtom();
            Atom oL = m.terminationAtom();
            iterator.reset(iF,oF,oL);
        }  
    }

    // Iterator for all intramolecular atom pairs within species
    public static class AllSpeciesIntraPairs implements AtomPairIterator.S {
        private boolean mLast;
        private Molecule m, lastM;
        private final AllMoleculeIntraPairs mPairs;
        public AllSpeciesIntraPairs(PhaseSpace ps) {  //constructor
            mPairs = new AllMoleculeIntraPairs(ps);  
            mLast = true;
        }
        public AllSpeciesIntraPairs(PhaseSpace ps, Species s) {  //constructor
            mPairs = new AllMoleculeIntraPairs(ps);  
            reset(s);
        }  
        public final boolean hasNext() {return (mPairs.hasNext() || !mLast);}
        public final AtomPair next() {
            if(!hasNext()) return null;
            if(!mPairs.hasNext()) {
                m = m.nextMolecule();
                mPairs.reset(m);
                mLast = (m == lastM);
            }
            return mPairs.next();                
        }
        public void reset(Species s) {
            if(s.getNAtomsPerMolecule() == 1) {m = null;} //no intramolecular atoms pairs, so jump to last molecule
            else {m = s.firstMolecule;}
            mPairs.reset(m);
            lastM = s.lastMolecule;
            mLast = (m==null);  //m null if only 1 atom/molecule, or no molecules in species;
        }
    }
    
    // Iterates over all intermolecular atom pairs between or within species
    public static class AllSpeciesInterPairs implements AtomPairIterator.S, AtomPairIterator.SS {
        private final AtomPairIterator.A fullIterator;  //make iterators in constructor and re-use them
        private final AtomPairIterator.A halfIterator;
        private AtomPairIterator.A currentIterator;
        private AtomPair thisPair, nextPair;
        private boolean hasNext, checkIntra;
        public AllSpeciesInterPairs(PhaseSpace ps) {  //constructor
            halfIterator = ps.makePairIteratorHalf();
            fullIterator = ps.makePairIteratorFull();
            hasNext = false;
        }  
        public AllSpeciesInterPairs(PhaseSpace ps, Species s) {  //constructor
            halfIterator = ps.makePairIteratorHalf();
            fullIterator = ps.makePairIteratorFull();
            reset(s);}  
        public AllSpeciesInterPairs(PhaseSpace ps, Species s1, Species s2) { //constructor
            halfIterator = ps.makePairIteratorHalf();
            fullIterator = ps.makePairIteratorFull();
            reset(s1, s2);}  
        public final boolean hasNext() {return hasNext;}
        public final AtomPair next() {
            thisPair = nextPair;  //save for return value, then prepare for next time around
            if(!currentIterator.hasNext()) {hasNext = false;}  //all done
            else {
                nextPair = currentIterator.next();
                if(checkIntra) { //check that not on same molecule; keep getting pair until satisfied
                    while(nextPair.atom1().parentMolecule == nextPair.atom2().parentMolecule) {
                        if(!currentIterator.hasNext()) {hasNext=false; break;}
                        nextPair = currentIterator.next();
                    }
                }
            }
            return thisPair;
        }
        public void reset(Species s) {reset(s,s);}
        public void reset(Species s1, Species s2) {
            Atom oF, oL, iF, iL;
            if(s1==null || s2==null) {hasNext=false; return;}
            oF = s1.firstAtom();                  
            if(oF==null) {hasNext=false; return;}    //no atoms in species 1
            iL = s2.lastAtom();
            if(s1 == s2) {  //same species
                if(oF==iL) {hasNext=false; return;}  //only one atom in species
                iF = oF.nextMoleculeFirstAtom();
                oL = iL.previousMoleculeLastAtom();
                checkIntra = (s1.getNAtomsPerMolecule() > 1);  //don't have to check for intramolecular pair if only one atom/molecule
                currentIterator = halfIterator;
                currentIterator.reset(iL,oF,oL);
            }
            else {          //different species
                if(iL==null) {hasNext=false; return;}  //no atoms in species 2
                iF = s2.firstAtom();                   //won't be null if iL!=null
                oL = s1.lastAtom();
                checkIntra = false;                    //different species, different molecules
                currentIterator = fullIterator;
                currentIterator.reset(iF,iL,oF,oL);
            }
            nextPair = currentIterator.next();
            hasNext = true;
        }
    }

}