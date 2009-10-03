package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IVectorMutable;

public interface IMoleculePositioned extends IMolecule {

    /**
     * Returns the position of the IMolecule.  Modifying the returned IVector
     * will alter the IMolecule's position, but not any positions of the
     * IMolecule's child IAtoms.
     */
    public IVectorMutable getPosition();
}
