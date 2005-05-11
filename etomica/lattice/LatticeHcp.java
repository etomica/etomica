package etomica.lattice;
import etomica.Default;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.space3d.Space3D;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

public class LatticeHcp extends LatticeCrystal {
    
    /**
     * Cubic hcp crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public LatticeHcp() {
        this(Default.ATOM_SIZE);
    }
    
    public LatticeHcp(double latticeConstant) {
        this(new PrimitiveHexagonal(Space3D.INSTANCE));
        primitive = (PrimitiveHexagonal)crystal.getLattice().getPrimitive();
        primitive.setA(latticeConstant);
        primitive.setC(Math.sqrt(8.0/3.0)*latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private LatticeHcp(PrimitiveHexagonal primitive) {
        super(new Crystal(primitive, new BasisHcp(primitive)));
    }
    
    /**
     * Returns the primitive the determines the lattice constant.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveHexagonal primitive() {
        return primitive;
    }

    public String toString() {return "Hcp";}
    
    private PrimitiveHexagonal primitive;
}//end of CrystalHcp