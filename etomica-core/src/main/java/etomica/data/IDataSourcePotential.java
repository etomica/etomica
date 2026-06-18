/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.data;

import etomica.potential.compute.PotentialCallback;

public interface IDataSourcePotential extends IDataSource {

    default boolean needsForces() {
        return false;
    }

    default boolean needsPairCallback() {
        return false;
    }

    default PotentialCallback getPotentialCallback() {
        return null;
    }

    void doCallComputeAll(boolean callComputeAll);
}
