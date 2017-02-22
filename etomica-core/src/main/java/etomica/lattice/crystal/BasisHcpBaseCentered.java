/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 2-atom basis for an hcp crystal 
 *   (base-centered or side-centered)
 *  
 *
 * @author Tai Boon Tan
 */
public class BasisHcpBaseCentered extends Basis {
    
    public BasisHcpBaseCentered() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
        new Vector3D(0.0, 0.0, 0.0),
        new Vector3D(0.5, 0.5, 0.0),
    };

    private static final long serialVersionUID = 1L;
}
