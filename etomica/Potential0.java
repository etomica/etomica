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
                    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential0Calculation) ) return;
        ((Potential0Calculation)pc).calculate(this); 
    }
        
}//end of Potential0



