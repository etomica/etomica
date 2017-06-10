/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.space.Space;

/**
 * 3D Hexagonal primitive with a 2-site hcp basis.
 */
public class LatticeHcp extends BravaisLatticeCrystal {
    
    /**
     * Hexagonal hcp crystal with a lattice constant that gives a
     * maximum-density structure for spheres of unit. 
     * <p>
     * Use scaleBy method if desired to make lattice constant give
     * maximum density for another sphere size.
     */
    public LatticeHcp(Space space) {
        this(space, 1.0);
    }
    
    public LatticeHcp(Space space, double latticeConstant) {
        super(new PrimitiveHexagonal(space, latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant), new BasisHcp());
        if(space.D() != 3) {
            throw new IllegalArgumentException("LatticeCubicHcp requires a 3-D space");
        }
    }

    /**
     * Rescales the lattice by the given factor. Multiplies the lattice constants
     * by the given value.
     */
    public void scaleBy(double scaleFactor) {
        primitive.scaleSize(scaleFactor);
    }

    public String toString() {return "Hcp";}
}
