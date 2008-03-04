package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IVector;

public interface IPotentialTorque extends PotentialSoft {

    public IVector[][] gradientAndTorque(IAtomSet atoms);
}
