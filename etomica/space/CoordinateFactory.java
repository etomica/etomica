package etomica.space;


/**
 * Interface for classes that make an Atom's coordinate field.  This is
 * used by the AtomFactory to make the coordinate that is passed to the
 * Atom's constructor.  The factory may be configured to generate coordinates
 * that do or do  not include velocities.
 */
public interface CoordinateFactory {

    /**
     * Returns a new coordinate instance.
     */
    public ICoordinate makeCoordinate();
    
    /**
     * Sets a flag that indicates if the Coordinate made by this factory
     * has velocity components.  If set to true, Coordinates returned by
     * makeCoordinate will have vectors to represent the velocity.
     */
    public void setKinetic(boolean kinetic);

    /**
     * Returns the value of the flag that indicates if the Coordinate
     * made by this factory has velocity components.
     */
    public boolean isKinetic();
    
}
