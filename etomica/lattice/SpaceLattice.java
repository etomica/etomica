package etomica.lattice;

import etomica.space.Space;



/**
 * Marker interface indicating that sites are instances of Space.Vector
 */

public interface SpaceLattice extends AbstractLattice {

    public Space getSpace();
    
    public double[] getLatticeConstants();
}
