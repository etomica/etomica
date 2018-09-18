/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Created by alex on 5/29/17.
 */
public class Vector3DTest {

    private static final double DELTA = 0.000000001;

    @Test
    public void testMod() throws Exception {
        Vector3D v1 = new Vector3D(-1, 7.5, 3);
        Vector3D v2 = new Vector3D(5.1, 5.1, 5.1);

        v1.mod(v2);

        double[] v = new double[3];
        v1.assignTo(v);
        assertArrayEquals(new double[]{4.1, 2.4, 3}, v, DELTA);

        v1 = new Vector3D(-7.5, 20.3, 0);

        v1.mod(v2);

        assertEquals(2.7, v1.x, DELTA);
        assertEquals(5, v1.y, DELTA);
        assertEquals(0, v1.z, DELTA);
    }

}
