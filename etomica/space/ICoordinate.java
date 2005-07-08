package etomica.space;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public interface ICoordinate {
    
    /**
     * Sets all parameters of this coordinate equal to those
     * of the given coordinate.
     * 
     * @param coord a coordinate of the same type as this coordinate
     */
    public void E(ICoordinate coord);

    public Vector position();

}