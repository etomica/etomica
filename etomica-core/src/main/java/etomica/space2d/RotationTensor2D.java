/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import etomica.space.RotationTensor;

public class RotationTensor2D extends Tensor2D implements RotationTensor {

    private static final long serialVersionUID = 1L;

    public RotationTensor2D() {
        super();
        reset();
    }

    public void reset() {
        xx = 1.0; xy = 0.0;
        yx = 0.0; yy = 1.0;
    }
    
    //FIXME should have a setOrientation

    /**
     *
     * @param i ignored
     * @param theta angle of rotation, in radians
     */
    public void setAxial(int i, double theta) {
        double st = Math.sin(theta);
        double ct = Math.cos(theta);
        xx = ct; xy=-st;
        yx = st; yy=ct;
    }
}
