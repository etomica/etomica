/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;
import etomica.lattice.crystal.BasisCubicBcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.space.Space;

/**
 * Cubic primitive with a 2-site bcc basis.
 */

public class LatticeCubicBcc extends BravaisLatticeCrystal implements CubicLattice {
    
	/**
	 * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of unit size.
     * <p>
     * Use scaleBy method if desired to make lattice constant give
     * maximum density for another sphere size.
	 */
    public LatticeCubicBcc(Space space) {
        this(space, 2.0/Math.sqrt(3.0));
    }
    
	public LatticeCubicBcc(Space space, double latticeConstant) {
		this(new PrimitiveCubic(space, latticeConstant));
        if(space.D() != 3) {
            throw new IllegalArgumentException("LatticeCubicBcc requires a 3-D space");
        }
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
	 * new BasisCubicBcc (which needs the new primitive) to super.
	 */	
	private LatticeCubicBcc(PrimitiveCubic primitive) {
		super(primitive, new BasisCubicBcc());
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
     * Returns "Bcc".
     */
    public String toString() {return "Bcc";}

    private static final long serialVersionUID = 1L;
}
