package etomica.space;

import etomica.atom.IAtomPositioned;

/**
 * Interface for a coordinate that includes an orientation.
 */
public interface ICoordinateAngular extends IAtomPositioned {

    public Orientation getOrientation();
}