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
//        public void reset();
        
        public interface SS extends Iterator {public void reset(Species.Agent s1, Species.Agent s2);}
        public interface S extends Iterator {public void reset(Species.Agent s);}
        public interface MS extends Iterator {public void reset(Molecule m, Species.Agent s);}
        public interface MM extends Iterator {public void reset(Molecule m1, Molecule m2);}
        public interface M extends Iterator {public void reset(Molecule m);}
        public interface A extends Iterator {  //core iterator, found in each Phase subclass
            public void reset(Atom a1, Atom a2, Atom a3); 
            public void reset(Atom a1, Atom a2, Atom a3, Atom a4);
            public void allDone();
            public void reset(Atom a, boolean intra);
        }
        
        /**
         * Iterator for all atom pairs in a phase
         * Default is to do inter and intra pairs; this may be overridden using reset method to do
         * only intermolecular pairs
         */
         public static class P implements Iterator {
            private final Iterator.A apiUp;
            private final Atom.Iterator atomUp;
            private boolean intra;
            private boolean hasNext;
            private AtomPair thisPair, nextPair;
            public P(Phase p) {
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
         * Iterator for all intermolecular atom pairs in a phase
         */
  /*      public static class P implements Iterator {
            private final Iterator.AMAM iteratorAMAM;
            private final Phase phase;
            private boolean hasNext;
            private int index1;
            private Potential2 p2;
            private Species.Agent s1, s2;
            private AtomPair nextPair, thisPair;
            public P(Phase ps) {iteratorAMAM = new AMAM(ps); phase = ps; hasNext = false;} //constructor
            public boolean hasNext() {return hasNext;}
            public AtomPair next() {
                thisPair = nextPair;
                if(iteratorAMAM.hasNext()) {               //another atom pair for this pair of species
                    nextPair = iteratorAMAM.next();
                    nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);
                }
                else {                                     //no more atom pairs for these two species
                    s2 = s2.nextSpecies();                 //increment inner loop over species
                    if(s2 == null) {                       //inner loop at end
                        s1 = s1.nextSpecies();             //increment outer loop over species
                        if(s1 == null) {hasNext = false;}  // outer loop at end --- all done
                        else {                             // restart inner loop                           
                            s2 = s1;                       // inner loop begins at current outer loop value (don't double-count species interactions)
                            index1 = s1.parentSpecies().speciesIndex;
                            p2 = phase.parentSimulation.potential2[index1][s2.parentSpecies().speciesIndex];
                            hasNext = iteratorAMAM.hasNext();
                            if(hasNext) {
                                nextPair = iteratorAMAM.next();   //assumes next species has a pair (this needs work)
                                nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);  
                            }
                        }
                    }
                    else {                                 // inner loop not done
                        p2 = phase.parentSimulation.potential2[index1][s2.parentSpecies().speciesIndex];
                        iteratorAMAM.reset(s1,s2);
                        hasNext = iteratorAMAM.hasNext();
                        if(hasNext) {
                            nextPair = iteratorAMAM.next();   //assumes next species has a pair (this needs work)
                            nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);  
                        }
                    }
                }
                return thisPair;
            } //end of next()
            public void reset() {
                s1 = phase.firstSpecies();
                s2 = s1;
                index1 = s1.parentSpecies().speciesIndex;
                p2 = phase.parentSimulation.potential2[index1][s2.parentSpecies().speciesIndex];
                iteratorAMAM.reset(s1,s2);
                hasNext = iteratorAMAM.hasNext();
                if(hasNext) {
                    nextPair = iteratorAMAM.next();   //assumes next species has a pair (this needs work)
                    nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);  
                }
            }
        }
  */          
        
        /**
         * Iterator for all atoms in a molecule with all atoms in a phase
         * The molecule may or may not be in the phase
         * Intramolecular pairs are not generated
         */
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
            public boolean hasNext() {return hasNext;}
            public AtomPair next() {
                thisPair = nextPair;
                if(apiCurrent.hasNext()) {nextPair = apiCurrent.next();}
                else {
                    if(upDone) {hasNext = false;}  //all done
                    else {                         //switch to down iterator
                        apiCurrent = apiDown;
                        upDone = true;
                        hasNext = apiCurrent.hasNext();
                        if(hasNext) nextPair = apiCurrent.next();
                    }
                }
                return thisPair;
            }
            public void reset(Molecule m) {
                apiUp.reset(m.lastAtom(),false);
                apiDown.reset(m.firstAtom(),false);
                if(apiUp.hasNext()) {
                    apiCurrent = apiUp;
                    upDone = false;}
                else {
                    apiCurrent = apiDown;
                    upDone = true;
                }
                hasNext = apiCurrent.hasNext();
            }
        }
                    
 /*           
        public static class MP implements M {
            private final Iterator.FMAM iteratorFMAM;
            private final Phase phase;
            private Molecule molecule;
            private boolean hasNext;
            private int index;       
            private Potential2 p2;
            private Species.Agent s1;
            private AtomPair nextPair, thisPair;
            public MP(Phase ps) {iteratorFMAM = new FMAM(ps); phase = ps; hasNext = false;} //constructor
            public MP(Phase ps, Molecule m) {iteratorFMAM = new FMAM(ps,m); phase = ps; reset(m);}
            public boolean hasNext() {return hasNext;}
            public AtomPair next() {
                thisPair = nextPair;
                if(iteratorFMAM.hasNext()) {       //iterator has another for this species
                    nextPair = iteratorFMAM.next();
                    nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);
                }
                else {                            //species used up; increment
                    s1 = s1.nextSpecies();
                    if(s1 == null) {hasNext = false;}  // no more species
                    else {                             // next species
                        p2 = phase.parentSimulation.potential2[s1.parentSpecies().speciesIndex][index];
                        iteratorFMAM.reset(s1);
                        nextPair = iteratorFMAM.next();   //assumes next species has a pair (this needs work)
                        nextPair.potential = p2.getPotential(nextPair.atom1, nextPair.atom2);
                    }
                }
                return thisPair;
            }
            public void reset(Molecule m) {
                molecule = m;
                s1 = phase.firstSpecies();
                iteratorFMAM.reset(m,s1);
                index = m.parentSpecies.speciesIndex;
                p2 = phase.parentSimulation.potential2[s1.parentSpecies().speciesIndex][index];
                hasNext = iteratorFMAM.hasNext();
                nextPair = iteratorFMAM.next();    //assumes has a next
            }
            public void reset() {reset(molecule);}
        }
 */       
        /** Naming scheme for iterators:
         *    AM --> Pairs involve All Molecules in a species, which is passed to constructor and/or reset method
         *    FM --> Pairs involve a Fixed Molecule, which is passed to constructor and/or reset method
         *     Single designator (e.g., AM) indicates INTRA-molecule pairs
         *     Double designator (e.g., AMFM) indicates INTER-molecule pairs
         */

        /**
         * Iterator for all intramolecular atom pairs within species
         */
  /*      public static class AM implements S {   //AllSpeciesIntraPairs
            private boolean mLast;
            private Molecule m, lastM;
            private FM mPairs;
            private Species.Agent species;
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
                species = s;
                if(s.parentSpecies().getAtomsPerMolecule() == 1) {m = null;} //no intramolecular atoms pairs, so jump to last molecule
                else {m = s.firstMolecule;}
                mPairs.reset(m);
                lastM = s.lastMolecule;
                mLast = (m==null);  //m null if only 1 atom/molecule, or no molecules in species;
            } 
            public void reset() {reset(species);}
        }
        
        /**
         * Iterator for all (intramolecular) atom pairs in a single molecule
         */
 /*       public static class FM implements M {  //AllMoleculeIntraPairs
            private final Iterator.A iterator;
            private Molecule molecule;
            public FM(Phase ps) {iterator = ps.makePairIteratorHalf();} //constructor
            public FM(Phase ps, Molecule m) {  //constructor
                iterator = ps.makePairIteratorHalf(); 
                reset(m);
            } 
            public boolean hasNext() {return iterator.hasNext();}
            public AtomPair next() {return iterator.next();}
            public final void reset(Molecule m) {
                molecule = m;
                if(m == null || m.nAtoms < 2) {iterator.allDone(); return;}
                Atom oF = m.firstAtom();
                Atom iF = oF.nextAtom();
                Atom oL = m.terminationAtom();
                iterator.reset(iF,oF,oL);
            } 
            public final void reset() {reset(molecule);}
        }

        /**
         * Iterates over all intermolecular atom pairs between a fixed molecule and all molecules in a species
         * Handles both cases in which molecule is or is not a member of the species
         */
  /*      public static class FMAM implements MS, M, S {
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
            public void reset() {reset(molecule,species);}
        }
        
        /**
         *  Iterates over all intermolecular atom pairs between or within species
         */
 /*       public static class AMAM implements S, SS {  //AllSpeciesInterPairs
            private final Iterator.A fullIterator;  //make iterators in constructor and re-use them
            private final Iterator.A halfIterator;
            private Species.Agent species1, species2;
            private Iterator.A currentIterator;
            private AtomPair thisPair, nextPair;
            private boolean hasNext, checkIntra;
            public AMAM(Phase ps) {  //constructor
                halfIterator = ps.iterator.makePairIteratorHalf();
                fullIterator = ps.iterator.makePairIteratorFull();
                hasNext = false;
            }  
            public AMAM(Phase ps, Species.Agent s) {  //constructor
                halfIterator = ps.iterator.makePairIteratorHalf();
                fullIterator = ps.iterator.makePairIteratorFull();
                reset(s);}  
            public AMAM(Phase ps, Species.Agent s1, Species.Agent s2) { //constructor
                halfIterator = ps.iterator.makePairIteratorHalf();
                fullIterator = ps.iterator.makePairIteratorFull();
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
                species1 = s1; species2 = s2;
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
            public void reset() {reset(species1, species2);}
        }  //end of class AMAM
  */      
 // Iterates over pairs formed by given atom and all atoms from other molecules above it in list
 // If given atom is not in phase, it is considered the last atom, and no iterations are performed
        public static final class Up implements Iterator.A {
            private final AtomPair pair;
            private final Phase phase;
            private boolean hasNext;
            private Atom nextAtom;
            public Up(Phase p) {phase = p; pair = new AtomPair(p); hasNext = false;}
            public Up(Phase p, Atom a) {phase = p; pair = new AtomPair(p); reset(a,true);}
            public void reset(Atom a1, Atom a2, Atom a3, Atom a4) {}
            public void reset(Atom a1, Atom a2, Atom a3) {}
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
            public void reset(Atom a1, Atom a2, Atom a3, Atom a4) {}
            public void reset(Atom a1, Atom a2, Atom a3) {}
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
        
        //These iterators need some work to permit iLast = null
        
    //"Full" --> Each iteration of inner loop begins with same first atom
    static final class Full implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        Atom outer, inner;
        private Atom iFirst, iLast, oFirst, oLast;
        private boolean hasNext;
        public Full(Phase p) {  //null constructor
            pair = new AtomPair(p);
            hasNext = false;
        }  
        public Full(Phase p, Atom iF, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair(p);
            reset(iF,iL,oF,oL);
        }
        public void reset(Atom iL, Atom oF, Atom oL) {reset(oF,iL,oF,oL);}  //take inner and outer first atoms as same
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {
            if(iF == null || oF == null) {hasNext = false; return;}
            iFirst = iF; 
            iLast =  iL; 
            oFirst = oF; 
            oLast =  oL;
            inner = iFirst;
            outer = oFirst;
            hasNext = true;
        }
        public void reset(Atom a, boolean intra) {}  
        public AtomPair next() {
            if(!hasNext) {return null;}
            pair.atom1 = outer;   //atom1 should always be outer
            pair.atom2 = inner;
            pair.reset();
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextAtom(); inner = iFirst;} //advance outer, reset inner
            }
            else {inner = inner.nextAtom();}
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
        public void reset() {reset(iFirst, iLast, oFirst, oLast);}
    }
 /*       
    //"Half" --> Each iteration of inner loop begins with atom after (HalfUp) or before (HalfDown) outer loop atom
    static final class HalfUp implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        Atom outer, inner;
        private Atom iFirst, iLast, oFirst, oLast;
        private boolean hasNext;
        public HalfUp(Phase p) {
            pair = new AtomPair(p);
            hasNext = false;
        }
        public HalfUp(Phase p, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair(p);
            reset(iL,oF,oL);
        }
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
        public void reset(Atom iL, Atom oF, Atom oL) {
            if(oF == null) {hasNext = false; return;}
            iLast =  iL; 
            oFirst =  oF; 
            oLast =  oL;
            outer = oFirst;
            inner = outer.nextAtom;
            hasNext = (inner != null);
        }
        public AtomPair next() {
            pair.atom1 = outer;   //c1 should always be outer
            pair.atom2 = inner;
            pair.reset();
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.nextAtom; inner = outer.nextAtom; hasNext = (inner!=null);} //advance outer, reset inner
            }
            else {inner = inner.nextAtom;}
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
        public void reset() {reset(iLast, oFirst, oLast);}
    }
    
    
    static final class HalfDown implements simulate.AtomPair.Iterator.A {
        final AtomPair pair;
        Atom outer, inner;
        private Atom iFirst, iLast, oFirst, oLast;
        private boolean hasNext;
        public HalfDown(Phase p) {
            pair = new AtomPair(p);
            hasNext = false;
        }
        public HalfDown(Phase p, Atom iL, Atom oF, Atom oL) {  //constructor
            pair = new AtomPair(p);
            reset(iL,oF,oL);
        }
        public void reset(Atom a, boolean intra) {
        public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
        public void reset(Atom iL, Atom oF, Atom oL) {
            if(oF == null) {hasNext = false; return;}
            iLast =  iL; 
            oFirst =  oF; 
            oLast =  oL;
            outer = oFirst;
            inner = outer.previousAtom();
            hasNext = (inner != null);
        }
        public AtomPair next() {
            pair.atom1 = outer;   //c1 should always be outer
            pair.atom2 = inner;
            pair.reset();
            if(inner == iLast) {                                     //end of inner loop
                if(outer == oLast) {hasNext = false;}                //all done
                else {outer = outer.previousAtom(); inner = outer.previousAtom(); hasNext = (inner != null);} //advance outer, reset inner
            }
            else {inner = inner.nextAtom();}
            return pair;
        }
        public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
        public boolean hasNext() {return hasNext;}
        public void reset() {reset(iLast, oFirst, oLast);}
    }
 */       
    }  //end of interface Iterator
}  //end of  AtomPair