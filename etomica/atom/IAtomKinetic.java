package etomica.atom;

import etomica.api.IAtomPositioned;
import etomica.api.IVector;

/**
 * Interface for an atom that holds vectors for
 * position and velocity.
 */
public interface IAtomKinetic extends IAtomPositioned {

    public IVector getVelocity();
 
}
