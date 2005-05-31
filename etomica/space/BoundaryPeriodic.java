package etomica.space;

/**
 * Interface for a boundary class that can be set to be periodic in
 * different directions.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 31, 2005 by kofke
 */
public interface BoundaryPeriodic {

    /**
     * Returns a boolean that indicates whether the boundary is periodic
     * in the corresponding direction.  Normally the array returned here
     * is the internal representation of the periodicity, so changing the
     * values in it will affect the periodicity properties of the boundary.
     */
    public boolean[] getPeriodicity();
}