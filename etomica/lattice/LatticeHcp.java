package etomica.lattice;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.space3d.Space3D;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

public class LatticeHcp extends BravaisLatticeCrystal {
    
    /**
     * Cubic hcp crystal with a lattice constant that gives a
     * maximum-density structure for spheres of unit. 
     * <p>
     * Use scaleBy method if desired to make lattice constant give
     * maximum density for another sphere size.
     */
    public LatticeHcp() {
        this(1.0);
    }
    
    public LatticeHcp(double latticeConstant) {
        this(new PrimitiveHexagonal(Space3D.getInstance(), latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant));
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private LatticeHcp(PrimitiveHexagonal primitive) {
        super(primitive, new BasisHcp());
    }
    
    /**
     * Rescales the lattice by the given factor. Multiplies the lattice constants
     * by the given value.
     */
    public void scaleBy(double scaleFactor) {
        primitive.scaleSize(scaleFactor);
    }

    public String toString() {return "Hcp";}

    private static final long serialVersionUID = 1L;
}