package etomica;

/**
 * Loops through all molecule pairs given within a phase.
 *
 * @author andrew
 */
 
/* History 
 * 08/05/04 Created
 */

public final class ApiIntraspeciesAA extends ApiIntralistAA implements AtomPairIterator, AtomPairPhaseIterator, AtomPairSpeciesIterator {
    
	private Species species;

	public ApiIntraspeciesAA(Simulation sim, Species species) {
		this(sim, species, species);
    }
	
    public ApiIntraspeciesAA(Simulation sim, Species species1, Species species2) {
    	super(sim);
    	if (species1 != species2) throw new RuntimeException("species passed to ApiIntraspeciesAA must be the same");
    	this.species = species1;
    }
    
	public void all(Phase phase, IteratorDirective dummy, AtomPairActive action) {
		all(((AtomTreeNodeGroup)phase.getAgent(species).node).childList, dummy, action);
	}
	
	public void setBasis(Phase phase) {
		setBasis(phase.speciesMaster().atomList);
	}

	public void setSpecies(Species species) {
		setSpecies(species);
	}
	
	public void setSpecies(Species species1, Species species2) {
    	if (species1 != species2) throw new RuntimeException("species passed to ApiIntraspeciesAA.setSpecies must be the same");
		this.species = species1;
	}
	
}  //end of class ApiIntraspeciesAA
    
