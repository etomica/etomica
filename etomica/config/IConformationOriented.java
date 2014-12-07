package etomica.config;

import etomica.api.IAtomList;
import etomica.space.IOrientation;

public interface IConformationOriented extends IConformation {

    /**
     * Assign coordinates to the atoms in the given list using the
     * given orientation for the molecule.
     */
    public void initializePositions(IAtomList atomList, IOrientation orientation);

}
