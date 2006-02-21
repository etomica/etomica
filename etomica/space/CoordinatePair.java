package etomica.space;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;

/**
 * Holds a pair of Coordinate instances and defines methods for computing distances and 
 * other quantities of interest related to the pair.  A CoordinatePair instance is
 * held by classes that manipulate and perform calculations based on pairs of atoms
 * (e.g., forces, collision times, radial distribution functions, etc.).  
 * <p>
 * {@link etomica.space.CoordinatePairKinetic CoordinatePairKinetic} must be used to compute 
 * quantities related to the relative velocities of the pair.
 */
public class CoordinatePair implements java.io.Serializable {

    /**
     * Constructs instance suitable for operation on coordinates defined for the
     * given Space.
     */
    public CoordinatePair(Space space) {
        dr = space.makeVector();
    }

    /**
     * Sets the class that defines how nearest images are determined when
     * calculating distances.  The nearest image is the pair of images that
     * are closest when all periodic-boundary images are considered.
     * Most often the class given here will be an instance of Boundary.
     * If this method is not called, no nearest-image transformation
     * is performed by reset.
     */
    public void setNearestImageTransformer(NearestImageTransformer b) {
        this.nearestImageTransformer = b;
    }

    /**
     * Accessor method for the class that defines how nearest images are determined.
     */
    public NearestImageTransformer getNearestImageTransformer() {
        return nearestImageTransformer;
    }
    
    /**
     * Specifies the pair of atoms to which this CoordinatePair applies.  Invoking 
     * this method causes the nearest-image separation between the atoms to
     * be calculated.
     * 
     * @return the separation vector between the pair of atoms, with nearest-image
     * transformation taken into account.  Separation is defined as pair.atom1 - pair.atom0.
     */
    public Vector reset(AtomPair pair) {
        return reset(((AtomLeaf)pair.atom0).coord, ((AtomLeaf)pair.atom1).coord);
    }

    /**
     * Specifies the pair of coordinates to which this CoordinatePair applies.  Invoking 
     * this method causes the nearest-image separation between the coordinates to
     * be calculated.
     * 
     * @return the separation vector between the pair of coordinates, with nearest-image
     * transformation taken into account.  Separation is defined as coord2 - coord1.
     */
    public Vector reset(ICoordinate coord1, ICoordinate coord2) {
        c1 = (Coordinate)coord1;
        c2 = (Coordinate)coord2;
        return reset();
    }

    /**
     * Calculates and returns the separation vector between the most-recently specific
     * coordinate/atom pair.
     * 
     * @return the separation vector between the pair of coordinates, with nearest-image
     * transformation taken into account.  Separation is defined as coord2 - coord1, where
     * "2" and "1" are as given in the last call to reset with arguments.
     * @throws NullPointerException if no coordinates or atoms were previously specified
     * or if the nearestImageTransformer was not specified.
     */
    public Vector reset() {
        dr.Ev1Mv2(c2.r, c1.r);
        nearestImageTransformer.nearestImage(dr);
        return dr;
    }

    /**
     * Returns the squared distance between the previously specified coordinates (with 
     * nearest-imaging considered).  Uses the separation calculated at the most recent
     * call to reset.  Returns zero if no atoms/coordinates were previously specified.
     */
    public final double r2() {
        return dr.squared();
    }

    /**
     * Returns the separation vector between the previously specified coordinates (with 
     * nearest-imaging considered).  Uses the separation calculated at the most recent
     * call to reset.  Returns a zero vector if no atoms/coordinates were previously specified.
     */
   public final Vector dr() {
        return dr;
    }

   /**
    * Translates the previously-specified position vectors by the amount rDelta in a
    * direction away from each other along the line connecting them.  Typically used
    * to separate or bring together by a small amount a pair of atoms having just 
    * completed a hard collision (to ensure they're on the right side of the collision potential).
    */
    public void nudge(double rDelta) {
        c1.position().PEa1Tv1(-rDelta, dr);
        c2.position().PEa1Tv1(+rDelta, dr);
    }
    
    protected Coordinate c1;
    protected Coordinate c2;
    protected final Vector dr;
    protected NearestImageTransformer nearestImageTransformer;
    private static final long serialVersionUID = 1L;
}
