package etomica.api;


/**
 * Interface for an atom that holds vectors for velocity.
 */
public interface IAtomKinetic extends IAtom {

    public IVectorMutable getVelocity();
 
}
