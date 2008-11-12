package etomica.api;


/**
 * Interface for an atom that holds vectors for
 * position and velocity.
 */
public interface IAtomKinetic extends IAtomPositioned {

    public IVector getVelocity();
 
}
