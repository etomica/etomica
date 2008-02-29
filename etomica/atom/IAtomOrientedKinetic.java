package etomica.atom;

import etomica.api.IVector;

/**
 * Interface for an Atom that has a position, orientation, velocity and angular
 * velocity.
 */
public interface IAtomOrientedKinetic extends IAtomKinetic, IAtomOriented {

    //XXX angular velocity is not a vector.  enjoy!
    public IVector getAngularVelocity(); //angular velocity vector in space-fixed frame
    
}