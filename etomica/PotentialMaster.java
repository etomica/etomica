package etomica;

import etomica.PotentialGroup.PotentialLinker;

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
    
	private PotentialGroupLrc lrcMaster;

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
 		for (PotentialLinker link=first; link!= null; link=link.next) {
			((AtomsetIteratorDirectable)link.iterator).setDirective(id);
		}
		for (PotentialLinker link=firstGroup; link!=null; link=link.next) {
			((AtomsetIteratorDirectable)link.iterator).setDirective(id);
		}
        for(PotentialLinker link=first; link!=null; link=link.next) {
        	((AtomIteratorPhaseDependent)link.iterator).setPhase(phase);
        	pc.doCalculation(link.iterator, link.potential);
        }//end for
        for(PotentialLinker link=firstGroup; link!=null; link=link.next) {
        	((AtomIteratorPhaseDependent)link.iterator).setPhase(phase);
            ((PotentialGroup)link.potential).calculate(link.iterator, id, pc);
        }
    }//end calculate
    
    public void setSpecies(Potential potential, Species[] species) {
    	if (species.length == 0 || potential.nBody() < species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
    	AtomsetIterator iterator;
    	switch (potential.nBody()) {
    	    case 0: break;
    	    case 1:
    	    	iterator = new AtomIteratorMolecule();
    	    	iterator.setSpecies(species[0]);
    	    	break;
    	    case 2:
    	    	if(species.length == 1 || species[0] == species[1]) {
    	    		iterator = new ApiIntragroup();
    	    	} else {
    	    		iterator = new ApiIntergroup();
    	    	}
    	    	AtomIteratorMolecule aiOuter = new AtomIteratorMolecule();
    	    	aiOuter.setSpecies(species[0]);
    	    	AtomIteratorNeighbor aiInner = new AtomIteratorNeighbor();
    	    	aiInner.setSpecies(species[species.length-1]);
    	    	iterator = new ApiInnerVariable(aiOuter, aiInner);
    	    	break;
    	    default:
    	    	throw new IllegalArgumentException("Handling potentials with more than 2-body interactions not implemented");
    	}
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
	
}//end of PotentialMaster
    