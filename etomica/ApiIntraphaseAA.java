package etomica;

/**
 * Loops through all molecule pairs given within a phase.
 *
 * @author andrew
 */
 
/* History 
 * 08/05/04 Created
 */

public final class ApiIntraphaseAA extends ApiIntralistAA implements AtomPairIterator, AtomPairPhaseIterator {
    
    public ApiIntraphaseAA(Simulation sim) {
    	super(sim);
    }
    
    public ApiIntraphaseAA(Simulation sim, Phase phase) {
    	this(sim);
    	setBasis(phase);
    }
    
	public void all(Phase phase, IteratorDirective dummy, AtomPairActive action) {
		all(phase.speciesMaster().atomList, dummy, action);
	}
	
	public void setBasis(Phase phase) {
		setBasis(phase.speciesMaster().atomList);
	}

}  //end of class ApiIntraspeciesAA
    
