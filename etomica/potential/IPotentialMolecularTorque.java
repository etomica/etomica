package etomica.potential;

import etomica.api.IMoleculeList;
import etomica.api.IVector;

public interface IPotentialMolecularTorque extends PotentialMolecularSoft {

    public IVector[][] gradientAndTorque(IMoleculeList atoms);
}
