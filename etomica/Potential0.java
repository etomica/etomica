package etomica; 

/**
 * Potential that does not depend on any atom positions.
 * Typically used to implement long-range corrections for potential truncation.
 * Potential thus depends on phase parameters, such as the number of molecules and the volume.
 *
 * @author David Kofke
 */

public abstract class Potential0 extends Potential {
  
    public static String VERSION = "Potential0:01.07.26/"+Potential.VERSION;
    
    public Potential0(SimulationElement parent) {
        super(0, parent);
    }
    
    public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
    	calculate(((SpeciesMaster)basis).node.parentPhase(), id, pc);                
    }
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
		pc.set(this).actionPerformed(phase);
    }
    
    public abstract double energy(Phase phase);
    
    public void setSpecies(Species[] species) {
    	throw new etomica.exception.MethodNotImplementedException();
//        switch (species.length) {
//            case 1: ;
//                    break;
//            default: throw new IllegalArgumentException("Wrong number of species given in Potential0");
//        }
//        super.setSpecies(species);
    }
    
    /**
     * Sets the iterator for this potential, which in most cases does not
     * require any iterator.  Default is AtomSetIterator.NULL.
     * @see etomica.Potential#setIterator(AtomSetIterator)
     */
    public void setIterator(AtomSetIterator iterator) {
    	this.iterator = iterator;
    }
    
    /**
     * Returns the iterator last defined via the setIterator method.  Default is
     * AtomSetIterator.NULL if none was previously set.
     * @see etomica.Potential#getIterator()
     */
    public AtomSetIterator getIterator() {
    	return iterator;
    }

    private AtomSetIterator iterator = AtomSetIterator.NULL;
        
}//end of Potential0



