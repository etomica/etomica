package simulate;

public interface AtomPair {
    
    public void reset();
    public double r2();
    public double v2();
    public double vDotr();
    public void push(double impulse);  //impart equal and opposite impulse to momenta
    public void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
    public Atom atom1();
    public Atom atom2();
    
    /**
     * Interface for classes that generate atom pairs according to various criteria
     */
    public interface Iterator {
        
        public boolean hasNext();
        public AtomPair next();
        
        public interface SS extends Iterator {public void reset(Species.Agent s1, Species.Agent s2);}
        public interface S extends Iterator {public void reset(Species.Agent s);}
        public interface MS extends Iterator {public void reset(Molecule m, Species.Agent s);}
        public interface MM extends Iterator {public void reset(Molecule m1, Molecule m2);}
        public interface M extends Iterator {public void reset(Molecule m);}
        public interface A extends Iterator {  //core iterator, found in each Phase subclass
            public void reset(Atom a1, Atom a2, Atom a3); 
            public void reset(Atom a1, Atom a2, Atom a3, Atom a4);
            public void allDone();
        }
        
        /** Naming scheme for iterators:
         *    AM --> Pairs involve All Molecules in a species, which is passed to constructor and/or reset method
         *    FM --> Pairs involve a Fixed Molecule, which is passed to constructor and/or reset method
         *     Single designator (e.g., AM) indicates INTRA-molecule pairs
         *     Double designator (e.g., AMFM) indicates INTER-molecule pairs
         */

        /**
         * Iterator for all intramolecular atom pairs within species
         */
        public static class AM implements S {   //AllSpeciesIntraPairs
            private boolean mLast;
            private Molecule m, lastM;
            private FM mPairs;
            public AM(Phase ps) {  //constructor
                mPairs = new FM(ps);
                mLast = true;
            }
            public AM(Phase ps, Species.Agent s) {  //constructor
                mPairs = new FM(ps);
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
            public void reset(Species.Agent s) {
                if(s.parentSpecies().getAtomsPerMolecule() == 1) {m = null;} //no intramolecular atoms pairs, so jump to last molecule
                else {m = s.firstMolecule;}
                mPairs.reset(m);
                lastM = s.lastMolecule;
                mLast = (m==null);  //m null if only 1 atom/molecule, or no molecules in species;
            } 
        }
        
        /**
         * Iterator for all (intramolecular) atom pairs in a single molecule
         */
        public static class FM implements M {  //AllMoleculeIntraPairs
            private final Iterator.A iterator;
            public FM(Phase ps) {iterator = ps.makePairIteratorHalf();} //constructor
            public FM(Phase ps, Molecule m) {  //constructor
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

        /**
         * Iterates over all intermolecular atom pairs between a fixed molecule and all molecules in a species
         */
        public static class FMAM implements MS, M, S {
            private final Iterator.A iterator;
            private Molecule molecule;
            private Species.Agent species;
            private boolean hasNext;
            private AtomPair thisPair, nextPair;
            private Atom oF, oL, iF, iL;
            public FMAM(Phase ps) {
                iterator = ps.makePairIteratorFull();
            }
            public FMAM(Phase ps, Molecule m) {
                iterator = ps.makePairIteratorFull();
                reset(m);
            }
            public FMAM(Phase ps, Species.Agent s) {
                iterator = ps.makePairIteratorFull();
                reset(s);
            }
            public FMAM(Phase ps, Molecule m, Species.Agent s) {
                iterator = ps.makePairIteratorFull();
                reset(m,s);
            }
            public final boolean hasNext() {return hasNext;}            
            public final AtomPair next() {
                thisPair = nextPair;  //save for return value, then prepare for next time around
                if(iterator.hasNext()) {nextPair = iterator.next();}  //hasNext remains true
                else if(iL==species.lastAtom() || molecule==species.lastMolecule()) {hasNext = false;}  //all done
                else {  //m in s, begins loop over atoms following s
                    iF = molecule.lastAtom.nextAtom();  //first atom on molecule after m
                    iL = species.lastAtom();
                    iterator.reset(iF,iL,oF,oL);
                    hasNext = iterator.hasNext();
                    if(hasNext) {nextPair = iterator.next();}
                }
                return thisPair;
            }
            
            public void reset(Molecule m, Species.Agent s) {
                molecule = m; species = s;
                if(m==null || s==null || s.nMolecules==0) {hasNext = false; return;}
                oF = m.firstAtom;   //molecule atoms are outer loop
                oL = m.lastAtom;
                if(m.parentPhase == s.parentPhase && m.parentSpecies == s.parentSpecies()) {  //m is a molecule of s
                    if(s.nMolecules == 0) {hasNext = false; return;}  //m is the only molecule of s
                    if(m == s.firstMolecule()) {   //m is first molecule of s; start loop after m
                        iF = oL.nextAtom();
                        iL = s.lastAtom();
                    }
                    else {                         //m is not the first molecule of s; start loop before m
                        iF = s.firstAtom();
                        iL = oF.previousAtom();
                    }
                }
                else {                             //m is not in s; loop over all atoms in s
                    iF = s.firstAtom();
                    iL = s.lastAtom(); 
                }
                iterator.reset(iF,iL,oF,oL);
                hasNext = iterator.hasNext();
                if(hasNext) {nextPair = iterator.next();}
            }
            public void reset(Molecule m) {reset(m,species);}
            public void reset(Species.Agent s) {reset(molecule,s);}
        }
        
        /**
         *  Iterates over all intermolecular atom pairs between or within species
         */
        public static class AMAM implements S, SS {  //AllSpeciesInterPairs
            private final Iterator.A fullIterator;  //make iterators in constructor and re-use them
            private final Iterator.A halfIterator;
            private Iterator.A currentIterator;
            private AtomPair thisPair, nextPair;
            private boolean hasNext, checkIntra;
            public AMAM(Phase ps) {  //constructor
                halfIterator = ps.makePairIteratorHalf();
                fullIterator = ps.makePairIteratorFull();
                hasNext = false;
            }  
            public AMAM(Phase ps, Species.Agent s) {  //constructor
                halfIterator = ps.makePairIteratorHalf();
                fullIterator = ps.makePairIteratorFull();
                reset(s);}  
            public AMAM(Phase ps, Species.Agent s1, Species.Agent s2) { //constructor
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
            public void reset(Species.Agent s) {reset(s,s);}
            public void reset(Species.Agent s1, Species.Agent s2) {
                Atom oF, oL, iF, iL;
                if(s1==null || s2==null) {hasNext=false; return;}
                oF = s1.firstAtom();                  
                if(oF==null) {hasNext=false; return;}    //no atoms in species 1
                iL = s2.lastAtom();
                if(s1 == s2) {  //same species
                    if(oF==iL) {hasNext=false; return;}  //only one atom in species
                    iF = oF.nextMoleculeFirstAtom();
                    oL = iL.previousMoleculeLastAtom();
                    checkIntra = (s1.parentSpecies().getAtomsPerMolecule() > 1);  //don't have to check for intramolecular pair if only one atom/molecule
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
            } //end of reset(s1,s2)
        }  //end of class AMAM
    }  //end of interface Iterator
}  //end of interface AtomPair