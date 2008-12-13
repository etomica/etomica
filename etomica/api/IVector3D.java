package etomica.api;


/**
 * Interface for a 3D vector, which has a XE (cross-product) method
 *
 * @author Andrew Schultz
 */
public interface IVector3D extends IVector {

    /**
     * Sets the vector components equal to the given values.
     */
    public void E(double a, double b, double c);
}