package etomica; 

/**
 * Potential that does not depend on any atom positions.
 * Typically used to implement long-range corrections for potential truncation.
 * Potential thus depends on phase parameters, such as the number of molecules and the volume.
 * Phase for application of potential is identified with the set(Phase) or set(SpeciesMaster)
 * methods.
 *
 * @author David Kofke
 */

public abstract class Potential0 extends Potential {
  
    public static String VERSION = "Potential0:01.07.26/"+Potential.VERSION;
    
    public Potential0(PotentialGroup parent) {
        super(parent);
    }
                    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential0Calculation) ) return;
        ((Potential0Calculation)pc).calculate(this); 
    }
    
    //we do not define energy method to take a Phase argument.  Instead
    //phase is specified with set method.  This is done to minimize the
    //times phase is specified, since in many cases it never changes during simulation
    public abstract double energy();
    
    public abstract Potential set(Phase phase);
    
    public Potential set(Atom[] atoms) {
        if(atoms.length != 0) throw new IllegalArgumentException("Too many atoms in Potential0");
        return this;
    }
    public Potential set(SpeciesMaster speciesMaster) {
        return set(speciesMaster.node.parentPhase());
    }
    
    public void setSpecies(Species s) {
        species = s;
    }
    
    public void setSpecies(Species[] species) {
        switch (species.length) {
            case 1: setSpecies(species[0]);
                    break;
            default: throw new IllegalArgumentException("Wrong number of species given in Potential2");
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



