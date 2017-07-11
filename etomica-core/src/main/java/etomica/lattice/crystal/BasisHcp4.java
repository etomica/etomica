/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 4-atom basis for an hcp crystal.
 *
 * @author Andrew Schultz
 */
public class BasisHcp4 extends Basis {

    public BasisHcp4() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
        new Vector3D(0.0, 0.0, 0.0),
        new Vector3D(0.5, 0.5,0),
        new Vector3D(0.5,1.0/3.0,0.5),
        new Vector3D(0,5.0/6.0,0.5)
    };
}
