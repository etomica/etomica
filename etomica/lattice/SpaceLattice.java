package etomica.lattice;

import etomica.space.Space;



/**
 * Marker interface indicating that sites are instances of Space.Vector
 */

/*
 * History
 * Created on Jan 3, 2005 by kofke
 */
public interface SpaceLattice extends AbstractLattice {

    public Space getSpace();
}
