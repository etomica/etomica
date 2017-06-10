/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;

/**
 * General basis class that hold scaled coordinates of atoms within a unit cell.
 *
 * @author David Kofke
 */
public class Basis implements java.io.Serializable {
    

    /**
     * @param scaledCoordinates basis coordinates for the case in which the
     * primitive lattice constant (getSize) is unity.  Given instance is kept
     * for interal representation of basis, so changes to scaledCoordinates
     * will affect the basis.
     */
    public Basis(Vector[] scaledCoordinates) {
        this.scaledCoordinates = scaledCoordinates;
    }
    
    /**
     * Returns scaled coordinates
     */
    public Vector[] getScaledCoordinates() {
        return scaledCoordinates;
    }
    
    private static final long serialVersionUID = 1L;
    private final Vector[] scaledCoordinates;
    
}
