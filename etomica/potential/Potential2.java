package etomica.potential;

import etomica.space.ISpace;

/**
 * Potential acting on 2 atoms or molecules.
 * 
 * @author David Kofke
 */

public abstract class Potential2 extends Potential {

    /**
     * Constructs potential with given space.
     */
    public Potential2(ISpace space) {
        super(2, space);
    }

}
