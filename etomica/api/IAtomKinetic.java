package etomica.api;


/**
 * Interface for an atom that holds vectors for velocity.
 */
public interface IAtomKinetic extends IAtom {

    /**
     * Returns the velocity of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's velocity.
     */
    public IVectorMutable getVelocity();
 
}
