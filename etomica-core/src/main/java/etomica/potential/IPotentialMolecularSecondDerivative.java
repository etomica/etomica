/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.molecule.IMoleculeList;
import etomica.space.Tensor;

public interface IPotentialMolecularSecondDerivative extends
		IPotentialMolecularTorque {
	 public Tensor [] secondDerivative(IMoleculeList molecules);
}
