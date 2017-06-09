/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

/**
 * OrientationCalc interface that handles quaternions (as a a double array)
 * instead of an Orientation object.
 */
public interface OrientationCalcQuaternion {
    public void calcOrientation(IMolecule molecule, double[] quat);
    public void setOrientation(IMolecule molecule, double[] quat);
}
