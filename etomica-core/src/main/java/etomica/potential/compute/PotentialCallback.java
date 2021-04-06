/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.potential.Potential2Soft;
import etomica.space.Tensor;
import etomica.space.Vector;

public interface PotentialCallback {

    default boolean wantsHessian() {return false;}

    default void pairCompute(int i, int j, Vector dr, double[] u012) {}

    /**
     *
     * @param drij
     * @param fij  force on j due to i (equal to -force on i due to j)
     * @param tij  torque on j due to i
     * @param tji  torque on i due to j
     */
    default void pairComputeGeneral(Potential2Soft pij, IAtom atom1, IAtom atom2, Vector drij, Vector fij, Vector tij, Vector tji) {}

    default void pairComputeHessian(int i, int j, Tensor phi) {}
}
