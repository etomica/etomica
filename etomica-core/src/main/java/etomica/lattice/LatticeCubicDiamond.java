/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;
import etomica.lattice.crystal.BasisCubicDiamond;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.space.Space;

/**
 * Cubic primitive with a 4-site fcc basis, on which each site 
 * is a 2-site diamond basis.
 */

 /* History
  * 09/26/02 (DAK) new
  * 01/20/04 (DAK) revised constructors; added one taking atomFactory argument
  */
public class LatticeCubicDiamond extends BravaisLatticeCrystal implements CubicLattice {
    
    /**
     * Cubic bcc crystal with a lattice constant that gives a maximum-density
     * structure for spheres of unit size.
     * <p>
     * Use scaleBy method if desired to make lattice constant give maximum
     * density for another sphere size.
     */
    public LatticeCubicDiamond(Space space) {
        this(space, 4.0/Math.sqrt(3.0));
    }
    
    public LatticeCubicDiamond(Space space, double latticeConstant) {
        this(new PrimitiveCubic(space));
        if(space.D() != 3) {
            throw new IllegalArgumentException("LatticeCubicDiamond requires a 3-D space");
        }
        ((PrimitiveCubic)primitive).setSizeABC(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private LatticeCubicDiamond(PrimitiveCubic primitive) {
        super(primitive, new BasisCubicDiamond());
    }
    
    /**
     * The lattice constant is the size of the cubic primitive vectors
     * of the lattice underlying this crystal.
     */
    public void setLatticeConstant(double latticeConstant) {
        ((PrimitiveCubic)primitive).setSizeABC(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return ((PrimitiveCubic)primitive).getSizeABC();
    }
    
    /**
     * Rescales the lattice by the given factor. Multiplies the lattice constant
     * by the given value.
     */
    public void scaleBy(double scaleFactor) {
        primitive.scaleSize(scaleFactor);
    }

    /**
     * Returns "Diamond".
     */
    public String toString() {return "Diamond";}
    
    private static final long serialVersionUID = 1L;
}
