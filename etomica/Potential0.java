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
    
    public Potential set(SpeciesMaster speciesMaster) {
        return set(speciesMaster.parentPhase());
    }
    
    public Potential set(Atom a) {//throw exception or redesign?
        return this;
    }
    public Potential set(Atom a1, Atom a2) {//throw exception or redesign?
        return this;
    }
        
}//end of Potential0



