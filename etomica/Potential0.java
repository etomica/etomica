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
    
    public Potential0(PotentialGroup parent) {
        super(parent);
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
        switch (species.length) {
            case 1: ;
                    break;
            default: throw new IllegalArgumentException("Wrong number of species given in Potential0");
        }
    }
    /**
     * Returns an array of length 1 with the species to which this potential applies.
     * Returns null if no species has been set.
     */
    public Species[] getSpecies() {
        if(species == null) return null;
        else return new Species[] {species};
    }

    private Species species;
        
}//end of Potential0



