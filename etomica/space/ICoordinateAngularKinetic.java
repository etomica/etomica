package etomica.space;

import etomica.atom.IAtomKinetic;

/**
 * Interface for a dynamical coordinate that includes an orientation in
 * addition to a spatial position.  
 */
public interface ICoordinateAngularKinetic extends IAtomKinetic, ICoordinateAngular {
    
    public IVector getAngularVelocity(); //angular velocity vector in space-fixed frame
    
}