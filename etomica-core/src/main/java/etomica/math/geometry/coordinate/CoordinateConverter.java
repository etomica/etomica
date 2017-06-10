/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry.coordinate;

import etomica.space.Vector;

public class CoordinateConverter {

    public static void toSpherical(Vector v, double[] result) {
        switch (v.getD()) {
            case 1:
                result[0] = v.getX(0);
                break;
            case 2:
                result[0] = Math.sqrt(v.squared());
                result[1] = Math.acos(v.getX(1) / result[0]); //theta
                break;
            case 3:
                result[0] = Math.sqrt(v.squared());
                result[1] = Math.acos(v.getX(2) / result[0]); //theta
                result[2] = Math.atan2(v.getX(2), v.getX(1)); //phi
                break;
            default:
                throw new IllegalArgumentException("I'm impressed by your ability to make a "+v.getD()+" dimension vector.  Please try 1-3 next time");
        }
    }
}
