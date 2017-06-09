package etomica.potential;

import etomica.atom.IMoleculeList;
import etomica.space.Tensor;

public interface IPotentialMolecularSecondDerivative extends
		IPotentialMolecularTorque {
	 public Tensor [] secondDerivative(IMoleculeList molecules);
}
