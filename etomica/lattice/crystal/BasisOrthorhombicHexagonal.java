/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space2d.Vector2D;

/**
 * A 2D hexagonal 2-atom basis
 *
 * @author Andrew Schultz
 */
public class BasisOrthorhombicHexagonal extends Basis {
    
    /**
     * Makes a 2D hexagonal 2-atom basis
     */
    public BasisOrthorhombicHexagonal() {
        super(scaledPositions);
    }
    
    private static final Vector2D[] scaledPositions = new Vector2D[] {
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 0.5)
    };
    
    private static final long serialVersionUID = 1L;
}