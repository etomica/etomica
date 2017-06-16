package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Tensor;

public interface IPotentialAtomicSecondDerivative extends
		IPotentialTorque {
	 Tensor [] secondDerivative(IAtomList atoms);
}
