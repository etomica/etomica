package etomica.space;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;

/**
 * Holds a pair of Coordinate instances and defines methods for computing
 * distances, relative velocities, and other quantities of interest related to
 * the pair. A CoordinatePair instance is held by classes that manipulate and
 * perform calculations based on pairs of atoms (e.g., forces, collision times,
 * radial distribution functions, etc.).
 * <p>
 * It is expected that classes instantitating this would be employed only if the
 * simulation is dynamic, meaning that the Coordinates a constructed with
 * velocity fields. Exceptions will be thrown if this class is used otherwise.
 * <p>
 * Note that the reset method (defined by the parent class) and the resetV
 * method are distinct, and must be invoked separately for position and velocity
 * calculations, respectively. A typical situation will invoke reset(AtomPair)
 * followed by resetV() (with no argument).
 * 
 */
public class CoordinatePairKinetic extends CoordinatePair {

    /**
     * Constructs instance suitable for operation on coordinates defined for the
     * given Space.
     */
    public CoordinatePairKinetic(Space space) {
        super(space);
        dv = space.makeVector();
    }

    /**
     * Specifies the pair of atoms to which this CoordinatePair applies and
     * calculates the relative velocity between them. No calculations are
     * performed related to the atoms' positions. The reset methods (as opposed
     * to resetV) do not perform any velocity-related calculations; the resetV
     * methods must be called explicitly if relative velocities are of interest.
     * 
     * @return the relative velocity of the pair of atoms, defined as pair.atom1 -
     *         pair.atom0
     * @throws ClassCastException
     *             if the atoms' coordinates do not implement ICoordinateKinetic
     */
    public void resetV(AtomPair pair) {
        resetV(((AtomLeaf)pair.atom0).coord, ((AtomLeaf)pair.atom1).coord);
    }

    /**
     * Specifies the pair of coordinates to which this CoordinatePair applies
     * and calculates the relative velocity between them. No calculations are
     * performed related to the coordinate positions. The reset methods (as
     * opposed to resetV) do not perform any velocity-related calculations; one
     * of the resetV methods must be called explicitly if relative velocities
     * are of interest.
     * 
     * @return the relative velocity of the pair of atoms, defined as coord2-coord1
     * @throws ClassCastException
     *             if the atoms' coordinates do not implement ICoordinateKinetic
     */
    public void resetV(ICoordinate coord1, ICoordinate coord2) {
        c1 = (Coordinate) coord1;
        c2 = (Coordinate) coord2;
        resetV();
    }

    /**
     * Calculates and returns the relative velocity for the pair of coordinates
     * most recently specified via the reset or resetV methods that take
     * arguments.
     * 
     * @return the relative velocity of the most recently specified pair of
     *         atoms, defined as described in the other resetV methods.
     * @throws ClassCastException
     *             if the atoms' coordinates do not implement ICoordinateKinetic
     * @throws NullPointerException
     *             if no atoms/coordinates were previously specified
     */
    public void resetV() {
        dv.Ev1Mv2(((ICoordinateKinetic) c2).velocity(), ((ICoordinateKinetic) c1).velocity());
    }

    /**
     * Returns the relative-velocity vector most recently calculated via resetV.
     */
    public final Vector dv() {
        return dv;
    }

    private final Vector dv;
    private static final long serialVersionUID = 1L;

}