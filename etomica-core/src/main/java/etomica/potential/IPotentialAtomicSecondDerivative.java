package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IMoleculeList;
import etomica.space.Tensor;

public interface IPotentialAtomicSecondDerivative extends
		IPotentialTorque {
	 public Tensor [] secondDerivative(IAtomList atoms);
}
