package etomica.lattice;

import etomica.Space;



/**
 * Marker interface indicating that sites are instances of Space.Vector
 */

/*
 * History
 * Created on Jan 3, 2005 by kofke
 */
public interface SpaceLattice extends AbstractLattice {

    public Space space();
}
