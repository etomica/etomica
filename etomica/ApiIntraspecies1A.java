package etomica;

/**
 * Loops through molecule pairs given within a single species, formed
 * from 1 molecule and All its neighbors.
 *
 * @author David Kofke
 */
public final class ApiIntraspecies1A implements AtomPairIterator {
    
    private Species species; 
    
    //this should be a neighbor iterator
    private final AtomIterator atomIterator = new AtomIteratorSequential();
    
    private Atom molecule;
    private final IteratorDirective localDirective = new IteratorDirective();
    private final AtomPair pair;
    
    public ApiIntraspecies1A() {
        atomIterator.setAsNeighbor(true);
        pair = new AtomPair(Simulation.instance.space);
    }
    
    public void setBasis(Atom a1, Atom a2) {
        if((a1 != a2) && !(a1 instanceof SpeciesAgent && a2 instanceof SpeciesAgent))
            throw new IllegalArgumentException("Improper basis given to ApiIntraSpecies1A");
        species = a1.node.parentSpecies();
        atomIterator.setBasis(a1);
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {return atomIterator.size();}        
    
    public boolean hasNext() {return atomIterator.hasNext();}
    
    public void reset(IteratorDirective id) {
        localDirective.set(id.direction());
        reset(id.atom1());
    }
    
    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        if(molecule == null) {
//            hasNext = false;
            return;
        }
        localDirective.set(IteratorDirective.BOTH);
        atomIterator.reset(localDirective.set(molecule));
        pair.atom1 = molecule;
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator directive (default UP is
     * if none previously specified.
     */
    public void reset(Atom atom) {
        if(atom == null) {
//            hasNext = false;
            return;
        }
        molecule = atom.node.parentMolecule();
        atomIterator.reset(localDirective.set(molecule));
        pair.atom1 = molecule;
    }
    
    public AtomPair next() {
        pair.atom2 = atomIterator.next();
        pair.reset();
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairAction act) {
        
    }
    
}  //end of class AtomPairIterator
    
