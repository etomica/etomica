package etomica;

/**
 * Iterator for the molecules in a phase.  Loops over those atoms
 * that lie just below the species agents in the atom tree hierarchy.
 * Can be set to loop over all molecules in a phase, or only those 
 * for a given species; use setBasis(Phase) and setBasis(SpeciesAgent),
 * respectively.
 *
 * @author David Kofke
 * @since 02.02.16
 */

public class AtomIteratorMolecule implements AtomIterator {
    
    public AtomIteratorMolecule() {
        treeIterator.setIterationDepth(2);
        iterator = treeIterator;
    }
    
    /**
     * Returns a new iterator ready to iterate over the molecules of
     * the given phase.
     */
    public AtomIteratorMolecule(Phase phase) {
        this();
        setBasis(phase);
        reset();
    }
    
    /**
     * Returns a new iterator ready to iterate over the molecules of
     * the given species in the given phase.
     */
    public AtomIteratorMolecule(Phase phase, Species species) {
        this();
        setBasis(phase, species);
        reset();
    }

    /**
	 * Invokes all(Atom, IteratorDirective, AtomActive) method of this
	 * class, using given arguments if they are instances of the appropriate
	 * classes. Otherwise returns without throwing any exception.
	 * @see etomica.AtomSetIterator#all(AtomSet, IteratorDirective, AtomSetActive)
	 */
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}

	public void all(Phase phase, IteratorDirective dummy, final AtomActive action) {
		all(phase.speciesMaster, dummy, action);    
	}
	public void all(Atom basis, IteratorDirective dummy, final AtomActive action) {
		if(basis instanceof SpeciesMaster) treeIterator.all(basis, dummy, action);
		else if(basis instanceof SpeciesAgent) listIterator.all(basis, dummy, action);
		else throw new IllegalArgumentException("Error: AtomIteratorMolecule invoked with illegal basis (not SpeciesMaster or SpeciesAgent");
	}	   
    /**
     * Puts iterator in a state in which hasNext is false.
     */
    public void unset() {iterator.unset();}
    
    public boolean hasNext() {return iterator.hasNext();}
    
    public boolean contains(Atom atom) {return iterator.contains(atom);}
    
    public Atom reset(IteratorDirective id) {
        throw new RuntimeException("method reset(IteratorDirective) not implemented in AtomIteratorMolecule");
    }
    
    public Atom reset() {return iterator.reset();}
    
    /**
     * Returns the next atom in the iteration sequence.
     */
    public Atom next() {return iterator.next();}
    
    public void allAtoms(AtomAction act) {iterator.allAtoms(act);}
    
    /**
     * Defines basis for iteration.  Must be an instance of SpeciesAgent or SpeciesMaster.
     */
    public void setBasis(Atom atom) {
        if(atom instanceof SpeciesMaster) setBasis((SpeciesMaster)atom);
        else if(atom instanceof SpeciesAgent) setBasis((SpeciesAgent)atom);
        else throw new IllegalArgumentException("Attempt to set inappropriate basis for AtomIteratorMolecule");
    }
    
    /**
     * Sets for iteration over all molecules in phase of given species master.
     */
    public void setBasis(SpeciesMaster speciesMaster) {
        treeIterator.setRoot(speciesMaster);
        iterator = treeIterator;
    }
    
    /**
     * Sets up for iteration over molecules under given species agent.
     */
    public void setBasis(SpeciesAgent speciesAgent) {
        listIterator.setList(((AtomTreeNodeGroup)speciesAgent.node).childList);
        iterator = listIterator;
    }
    
    /**
     * Sets for iteration over all molecules in given phase.
     */
    public void setBasis(Phase phase) {setBasis(phase.speciesMaster);}
    
    /**
     * Sets for iteration over molecules of given species in given phase.
     */
    public void setBasis(Phase phase, Species species) {setBasis(species.getAgent(phase));}
        
    
    public Atom getBasis() {return iterator.getBasis();}
    
    public int size() {return iterator.size();}
    
    public void setAsNeighbor(boolean b) {
        throw new RuntimeException("setAsNeighbor not implemented in AtomIteratorMolecule");
    }
    
    private final AtomIteratorTree treeIterator = new AtomIteratorTree();
    private final AtomIteratorList listIterator = new AtomIteratorList();
    private AtomIterator iterator;
    
    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {
        
        Simulation sim = new Simulation();
        Simulation.instance = sim;
        Species species2 = new SpeciesSpheresMono();
        Species species1 = new SpeciesSpheres(3,3);
        Species species0 = new SpeciesSpheres(3,2);
        species0.setNMolecules(5);
        species1.setNMolecules(2);
        species2.setNMolecules(3);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
        
        AtomIteratorMolecule iterator = new AtomIteratorMolecule();
        System.out.println(iterator.hasNext());
        
        iterator.setBasis(phase);
        System.out.println(iterator.hasNext());
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        System.out.println();
        
        iterator.setBasis(phase, species2);
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        System.out.println();
    }//end of main
    
    
}//end of AtomIteratorMolecule
