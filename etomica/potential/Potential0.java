package etomica.potential; 

import etomica.space.Space;

/**
 * Potential that does not depend on any atom positions.
 * Typically used to implement long-range corrections for potential truncation.
 * Potential thus depends on phase parameters, such as the number of molecules 
 * and the volume.
 *
 * @author David Kofke
 */

public abstract class Potential0 extends Potential {
      
    public Potential0(Space space) {
        super(0, space);
    }
    
    /**
     * Returns zero.
     */
    public double getRange() {
        return 0.0;
    }
                        
}//end of Potential0



