package etomica;

/**
 * Iterator for the molecules in a phase.  Loops over those atoms
 * that lie just below the species agents in the atom tree hierarchy.
 * Can be set to loop over all molecules in a phase, or only those 
 * for a given species.
 *
 * @author David Kofke
 * @since 02.02.16
 */

public final class AtomIteratorMolecule extends AtomIteratorAdapter implements AtomsetIteratorSpeciesDependent {
    
    public AtomIteratorMolecule() {
        super(new AtomIteratorTree(2));
        treeIterator = (AtomIteratorTree)iterator;
        setPhase(null);
    }
    
    /**
     * Returns a new iterator ready to iterate over the molecules of
     * the given phase.
     */
    public AtomIteratorMolecule(Phase phase) {
        this();
        setPhase(phase);
    }
    
    /**
     * Returns a new iterator ready to iterate over the molecules of
     * the given species in the given phase.
     */
    public AtomIteratorMolecule(Phase phase, Species species) {
        this();
        set(phase, species);
    }
    
    public void set(Phase phase, Species species) {
    	this.species = species;
    	setPhase(phase);
    }
    
    public void setPhase(Phase phase) {
    	if(phase == null) treeIterator.setRoot(null);
    	else if(species == null) {
    		treeIterator.setIterationDepth(2);
    		treeIterator.setRoot(phase.speciesMaster);
    	} else {
    		treeIterator.setIterationDepth(1);
    		treeIterator.setRoot(phase.getAgent(species));
    	}
    }
    
    public void setSpecies(Species[] species) {
    	this.species = species[0];
    	setPhase(phase);
    }
    
    private final AtomIteratorTree treeIterator;
    private Species species;
    private Phase phase;
    
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
        
        iterator.setPhase(phase);
        System.out.println(iterator.hasNext());
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        System.out.println();
        
        iterator.set(phase, species2);
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        System.out.println();
    }//end of main
    
    
}//end of AtomIteratorMolecule
