package etomica.space;

/**
 * Interface for a Coordinate that includes a Vector that defines the atom's position.
 */
public interface ICoordinate {
    
    /**
     * Sets all parameters of this coordinate equal to those
     * of the given coordinate.
     * 
     * @param coord a coordinate of the same type as this coordinate
     */
    public void E(ICoordinate coord);

    /**
     * Returns the position vector.
     */
    public Vector getPosition();

}