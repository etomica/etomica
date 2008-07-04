package etomica.atom;

import etomica.api.IMolecule;

/**
 * OrientationCalc interface that handles quaternions (as a a double array)
 * instead of an Orientation object.
 */
public interface OrientationCalcQuaternion {
    public void calcOrientation(IMolecule molecule, double[] quat);
    public void setOrientation(IMolecule molecule, double[] quat);
}