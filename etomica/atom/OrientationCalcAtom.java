package etomica.atom;

import etomica.api.IMolecule;
import etomica.space3d.IOrientationFull3D;

/**
 * OrientationCalc implementation that handles a monotomic oriented molecule.
 *
 * @author Andrew Schultz
 */
public class OrientationCalcAtom implements OrientationCalc {

    public void calcOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        orientation.E(((IAtomOriented)molecule).getOrientation());
    }

    public void setOrientation(IMolecule molecule,
            IOrientationFull3D orientation) {
        ((IAtomOriented)molecule).getOrientation().E(orientation);
    }

}
