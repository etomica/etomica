package etomica;

/**
 * Loops through molecule pairs given within a single group, formed
 * from 1 molecule and All its neighbors.  Molecule is specified itself
 * or via one of its descendant atoms, via an iterator directive.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/04/02 (DAK) Changed reset(atom) to use isDescendedFrom rather than just
  * checking parent. Changes performed to correct problem uncovered in
  * "surfactant" module simulation
  * 
  * 01/02/03 (DAK) reset checks if basis is null
  * 01/27/03 (DAK) added "all" method, converted to class as part of redesign of
  * Potential
  */
 
public final class ApiIntragroup1A extends AtomPairIterator {
    
    public ApiIntragroup1A(Simulation sim) {
        pair = new AtomPair(sim.space);
        atomIterator = sim.iteratorFactory.makeIntragroupNbrIterator();
    }
    
	public void all(Atom basis, IteratorDirective id, AtomPairActive action) {
		Atom atom = id.atom1();
		if(atom == null || basis == null || action == null) return;
		AtomTreeNode referenceNode = atom.node.childWhereDescendedFrom(basis.node);
		if(referenceNode == null) return;
		AtomPairActive.InnerWrapper wrapper = action.innerWrapper();
		wrapper.pair.atom1 = referenceNode.atom;
		atomIterator.all(basis, id.set(referenceNode.atom), wrapper);
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairActive action) {
		throw new IllegalArgumentException("Error: AtomPairIterator not defined for a pair basis");
	}

    public void setBasis(Atom a1, Atom a2) {
        if(a1 != a2)
            throw new IllegalArgumentException("Improper basis given to ApiIntragroup1A");
        atomIterator.setBasis(a1);
        group = a1;
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
     * Resets the iterator using the current reference atom, or unsets if no
     * atom was previously given.
     */
    public void reset() {
        if(referenceAtom == null || group == null) atomIterator.unset();
        else atomIterator.reset(localDirective.set(IteratorDirective.BOTH).set(referenceAtom));
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator direction (default is BOTH
     * if none previously specified).
     */
    public void reset(Atom atom) {
        if(atom == null || group == null) {//group could be null if potential is "turned off"
            atomIterator.unset();
            return;
        }
        AtomTreeNode referenceNode = atom.node.childWhereDescendedFrom(group.node);
        if(referenceNode == null) atomIterator.unset();
        else {
            referenceAtom = referenceNode.atom;
            atomIterator.reset(localDirective.set(referenceAtom));
        }
        pair.atom1 = referenceAtom;
    }
    
    public AtomPair next() {
        pair.atom2 = atomIterator.next();
        pair.reset();
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     * Must have called reset(Atom) before invoking.
     */
    public void allPairs(AtomPairAction act) {
    	throw new etomica.exception.MethodNotImplementedException();
//        if(referenceAtom == null) return;
//        atomIterator.allAtoms(act.wrapper);
    }
    
    private final AtomIterator atomIterator;
    
    private Atom group;
    private Atom referenceAtom;
    private final IteratorDirective localDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomPair pair;
        
    /**
     * Method to test and demonstrate use of class.
     */
    public static void main(String[] args) {
        
        Simulation sim = new Simulation();
        int nMolecules = 5;
        SpeciesSpheresMono species = new SpeciesSpheresMono();
        species.setNMolecules(nMolecules);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
//        AtomList atomList = phase.speciesMaster.atomList;
        AtomList atomList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
        
        AtomPairIterator iterator = new ApiIntragroup1A(sim);
        Atom first = phase.speciesMaster.firstSpecies().firstMolecule();
        Atom last = phase.speciesMaster.lastSpecies().lastMolecule();
        Atom middle = null;
        if(nMolecules > 2) {
            do middle = phase.randomMolecule();
            while (middle == first || middle == last);
        }
        
        iterator.setBasis(phase.speciesMaster.firstSpecies(), phase.speciesMaster.lastSpecies());
        IteratorDirective.testSuitePair(iterator, first, middle, last);
    }
    
    
}  //end of class ApiIntragroup1A
    
