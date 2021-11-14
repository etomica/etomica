package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Tensor;

public interface IPotentialAtomicSecondDerivative {
	 Tensor [] secondDerivative(IAtomList atoms);
}
