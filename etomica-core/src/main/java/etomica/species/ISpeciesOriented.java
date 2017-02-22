/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.api.ISpecies;
import etomica.api.IVector;

public interface ISpeciesOriented extends ISpecies {

    /**
     * Returns the principle components of the moment of inertia of the
     * molecule within the body-fixed frame.  Do NOT modify the returned moment
     * of inertia returned.
     */
    public IVector getMomentOfInertia();

    public double getMass();

}