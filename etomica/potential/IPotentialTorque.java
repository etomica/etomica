package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IVector;

public interface IPotentialTorque extends PotentialSoft {

    public IVector[][] gradientAndTorque(IAtomList atoms);
}
