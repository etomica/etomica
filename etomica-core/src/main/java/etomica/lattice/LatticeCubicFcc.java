/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.space.Space;

/**
 * Cubic primitive with a 4-site fcc basis.
 * 
 */
public class LatticeCubicFcc extends BravaisLatticeCrystal implements CubicLattice {

    /**
     * Cubic fcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of unit size. 
     * <p>
     * Use scaleBy method if desired to make lattice constant give
     * maximum density for another sphere size.
     */
    public LatticeCubicFcc(Space space) {
        this(space, Math.sqrt(2.0));
    }
    
    public LatticeCubicFcc(Space space, double latticeConstant) {
        this(new PrimitiveCubic(space, latticeConstant));
        if(space.D() != 3) {
            throw new IllegalArgumentException("LatticeCubicFcc requires a 3-D space");
        }
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicFcc (which needs the new primitive) to super.
     */ 
    private LatticeCubicFcc(PrimitiveCubic primitive) {
        super(primitive, new BasisCubicFcc());
    }
    
    /**
     * Returns a new PrimitiveFcc instance corresponding to the fcc lattice. 
     * Changes to the returned instance have no effect on the lattice.
     */
    public PrimitiveFcc getPrimitiveFcc() {
        PrimitiveFcc p = new PrimitiveFcc(getSpace());
        p.setCubicSize(getLatticeConstant()/Math.sqrt(2.0));
        return p;
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
     * Returns "Fcc".
     */
    public String toString() {return "Fcc";}
    
    private static final long serialVersionUID = 1L;
}
