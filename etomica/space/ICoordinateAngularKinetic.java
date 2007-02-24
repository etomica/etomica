package etomica.space;

/**
 * Interface for a dynamical coordinate that includes an orientation in
 * addition to a spatial position.  
 */
public interface ICoordinateAngularKinetic extends ICoordinateKinetic, ICoordinateAngular {
    
    public IVector getAngularVelocity(); //angular velocity vector in space-fixed frame
    
}