package etomica.potential;

import etomica.atom.AtomSet;
import etomica.api.IVector;

public interface IPotentialTorque extends PotentialSoft {

    public IVector[][] gradientAndTorque(AtomSet atoms);
}
