/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IRandom;
import etomica.api.IVectorMutable;

/**
 * Vector interface for vectors having random number methods 
 */
public interface IVectorRandom extends IVectorMutable {

    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public void setRandomSphere(IRandom random);

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public void setRandomCube(IRandom random);

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public void setRandomInSphere(IRandom random);
}