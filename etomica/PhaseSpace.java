package simulate;

public abstract class PhaseSpace extends Component {

    
    public abstract PhaseSpaceCoordinate makeCoordinate();
    public abstract AtomPairIterator makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator makePairIteratorHalf(Atom iL, Atom oF, Atom oL);
    public abstract AtomPairIterator makePairIteratorFull();
    public abstract AtomPairIterator makePairIteratorHalf();
    

    // Iterates over all intermolecular atom pairs between species (which may be the same species)
    public class AllSpeciesInterPairs implements AtomPairIterator {
        private final AtomPairIterator fullIterator = makePairIteratorFull();  //make iterators here and re-use them
        private final AtomPairIterator halfIterator = makePairIteratorHalf();
        private AtomPairIterator currentIterator;
        private AtomPair thisPair, nextPair;
        private boolean hasNext;
        public AllSpeciesInterPairs() {hasNext = false;}  //constructor
        public AllSpeciesInterPairs(Species s1, Species s2) {reset(s1, s2);}  //constructor
        public final boolean hasNext() {return hasNext;}
        public final SpaceAtomPair next() {
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