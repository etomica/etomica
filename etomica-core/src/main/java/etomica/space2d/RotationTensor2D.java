/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import etomica.space.RotationTensor;

public class RotationTensor2D extends Tensor2D implements RotationTensor {
    public RotationTensor2D() {super(); reset();}

    public void reset() {
        xx = 1.0; xy = 0.0;
        yx = 0.0; yy = 1.0;
    }

    public void setAxial(int i, double theta) {
        double st = Math.sin(theta);
        double ct = Math.cos(theta);
        switch(i) {
            case 2: xx = ct; xy=-st;
                    yx = st; yy=ct;
                    return;
            default: throw new IllegalArgumentException();
        }
    }
    
    //FIXME should have a setOrientation
    
    private static final long serialVersionUID = 1L;
}
