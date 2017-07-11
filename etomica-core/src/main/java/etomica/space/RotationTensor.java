/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

/**
 * A tensor intended for use to perform rotations in space.  Includes methods
 * to set its components in terms of rotation angles. Methods are formulated
 * to ensure that application to a vector does not change magnitude of vector.
 */
public interface RotationTensor extends Tensor {
    
    /**
     * Sets tensor to a condition of no rotation.
     */
    public void reset();

    public void setAxial(int i, double theta);
}