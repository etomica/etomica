package etomica;

/**
 * Loops through molecule pairs given between two species, formed
 * from 1 molecule and All its neighbors.  Molecule is specified itself
 * or in terms of one of its descendant atoms, via an iterator directive.
 *
 * @author David Kofke
 */
 
public final class ApiInterspecies1A implements AtomPairIterator {
    
    private SpeciesAgent speciesAgent1; 
    private SpeciesAgent speciesAgent2;
    
    //this should be a neighbor iterator
    private final AtomIterator atomIterator = new AtomIteratorChildren();
    
    private Atom molecule;
    private final IteratorDirective localDirective = new IteratorDirective();
    private final AtomPair pair;
    
    public ApiInterspecies1A(Simulation sim) {
        pair = new AtomPair(sim.space);
    }
    
    /**
     * Identifies the species agents that are the parents of the molecules
     * to be given by the iterator.  The arguments must be different instances
     * of SpeciesAgent.
     */
    public void setBasis(Atom a1, Atom a2) {
        if((a1 == a2) || !(a1 instanceof SpeciesAgent && a2 instanceof SpeciesAgent))
            throw new IllegalArgumentException("Improper basis given to ApiInterSpecies1A");
        speciesAgent1 = (SpeciesAgent)a1;
        speciesAgent2 = (SpeciesAgent)a2;
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator, based
     * on the most recent specification of the atom given to reset.
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
            return;
        }
        atomIterator.reset();
        pair.atom1 = molecule;
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom.
     */
    public void reset(Atom atom) {
        if(atom == null) return;
        molecule = atom.node.parentMolecule();
        SpeciesAgent agent = (SpeciesAgent)molecule.node.parentGroup();
        if(agent == speciesAgent1) {
            atomIterator.setBasis(speciesAgent2);
        } else if(agent == speciesAgent2) {
            atomIterator.setBasis(speciesAgent1);
        } else {
            atomIterator.setBasis(null);
            atomIterator.reset();
            return;
        }
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
     * Not yet implemented.
     */
    public void allPairs(AtomPairAction act) {
        throw new RuntimeException("Method allPairs not implemented in ApiInterspecies1A");
    }
    
}  //end of class AtomPairIterator
    
