package etomica;

/**
 * Master potential that oversees all other potentials in the Hamiltonian.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation 
 * on to the contained potentials.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/13/02 (DAK) added removePotential method
  * 01/27/03 (DAK) many revisions as part of redesign of Potential
  * 06/15/03 (DAK) in makeBasis, for 2-species case, added check for whether the
  * two species are the same; if so uses one of them as basis, rather than
  * making a pair from them.
  */
public final class PotentialMaster extends PotentialGroup {
    
    public PotentialMaster(Simulation sim) {
        super(-1, sim);
    } 
    
	/**
	 * Returns the potential group that oversees the long-range
	 * correction zero-body potentials.
	 */
	 public PotentialGroupLrc lrcMaster() {
		if(lrcMaster == null) lrcMaster = new PotentialGroupLrc(this);
		return lrcMaster;
	 }

	 public void calculate(AtomsetIterator iterator, IteratorDirective id, PotentialCalculation pc) {
	 	throw new RuntimeException("Method inappropriate for PotentialMaster class");
	 }
	 
	//should build on this to do more filtering of potentials based on directive
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
    	if(!enabled) return;
    	Atom[] targetAtoms = id.targetAtoms();
    	boolean phaseChanged = (phase != mostRecentPhase);
    	mostRecentPhase = phase;
    	for(PotentialLinker link=first; link!=null; link=link.next) {
			if(phaseChanged) ((AtomsetIteratorPhaseDependent)link.iterator).setPhase(phase);
			((AtomsetIteratorTargetDependent)link.iterator).setTarget(targetAtoms);
        	pc.doCalculation(link.iterator, link.potential);
        }//end for
        for(PotentialLinker link=firstGroup; link!=null; link=link.next) {
        	if(phaseChanged) ((AtomsetIteratorPhaseDependent)link.iterator).setPhase(phase);
			((AtomsetIteratorTargetDependent)link.iterator).setTarget(targetAtoms);
           ((PotentialGroup)link.potential).calculate(link.iterator, id, pc);
        }
    }//end calculate
    
    public void setSpecies(Potential potential, Species[] species) {
    	if (species.length == 0 || potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
    	AtomsetIteratorMolecule iterator = new AtomsetIteratorMolecule(species);
    	addPotential(potential, iterator);
    }
 
	/**
	 * Convenient reformulation of the calculate method, applicable if the
	 * potential calculation performs a sum.  The method returns the
	 * summable potential calculation object, so that the sum can be accessed
	 * in-line with the method call.
	 */
   public final PotentialCalculation.Summable calculate(Phase phase, IteratorDirective id, PotentialCalculation.Summable pa) {
	   this.calculate(phase.speciesMaster, id, (PotentialCalculation)pa);
	   return pa;
   }	    
   
	private PotentialGroupLrc lrcMaster;
	private Phase mostRecentPhase = null;

	
}//end of PotentialMaster
    