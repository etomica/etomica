package simulate;

public abstract class PhaseSpace extends Component {

    
    public abstract PhaseSpace.Coordinate makeCoordinate();
    public abstract AtomPairIterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator.A makePairIteratorFull();
    public abstract AtomPairIterator.A makePairIteratorHalf();
    
    
    abstract class Coordinate {
        public abstract void setNextCoordinate(PhaseSpace.Coordinate c);
        public abstract void clearPreviousCoordinate();
        
        public final void setNextAtom(Atom a) {
            if(a == null) {setNextCoordinate(null);}
            else {setNextCoordinate(a.coordinate);}
        }
        public final void clearPreviousAtom() {clearPreviousCoordinate();}
    }
    
    // Iterator for all intramolecular atom pairs within species
    public class AllSpeciesIntraPairs implements AtomPairIterator.S {
        private boolean mLast;
        private Molecule m, lastM;
        AllMoleculeIntraPairs mPairs = new AllMoleculeIntraPairs();
        public AllSpeciesIntraPairs() {mLast = true;}
        public AllSpeciesIntraPairs(Species s) {reset(s);}  //constructor
        public final boolean hasNext() {return (mPairs.hasNext() || !mLast);}
        public final SpaceAtomPair next() {
            if(!hasNext()) return null;
            if(!mPairs.hasNext()) {
                m = m.getNextMolecule();
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
    
    // Iterator for all (intramolecular) atom pairs in a single molecule
    public class AllMoleculeIntraPairs implements AtomPairIterator.M {
        private final AtomPairIterator.A iterator = makePairIteratorHalf();
        public AllMoleculeIntraPairs() {}
        public AllMoleculeIntraPairs(Molecule m) {reset(m);}
        public boolean hasNext() {return iterator.hasNext();}
        public AtomPair next() {return iterator.next();}
        public final void reset(Molecule m) {
            if(m == null || m.nAtoms < 2) {iterator.allDone(); return;}
            Atom oF = m.firstAtom();
            Atom iF = oF.getNextAtom();
            Atom oL = m.terminationAtom();
            iterator.reset(iF,oF,oL);
        }  
    }
    
    // Iterates over all intermolecular atom pairs between or within species
    public class AllSpeciesInterPairs implements AtomPairIterator.S, AtomPairIterator.SS {
        private final AtomPairIterator.A fullIterator = makePairIteratorFull();  //make iterators here and re-use them
        private final AtomPairIterator.A halfIterator = makePairIteratorHalf();
        private AtomPairIterator currentIterator;
        private AtomPair thisPair, nextPair;
        private boolean hasNext;
        public AllSpeciesInterPairs() {hasNext = false;}  //constructor
        public AllSpeciesInterPairs(Species s) {reset(s);}  //constructor
        public AllSpeciesInterPairs(Species s1, Species s2) {reset(s1, s2);}  //constructor
        public final boolean hasNext() {return hasNext;}
        public final AtomPair next() {
            thisPair = nextPair;  //save for return value, then prepare for next time around
            if(!currentIterator.hasNext()) {hasNext = false;}  //all done
            else {
                nextPair = currentIterator.next();
                if(checkIntra) { //check that not on same molecule; keep getting pair until satisfied
                    while(nextPair.atom1.parentMolecule == nextPair.atom2.parentMolecule) {
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