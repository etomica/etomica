package etomica.lattice.crystal;
import etomica.Default;
import etomica.Space3D;
import etomica.lattice.BravaisLattice;
import etomica.lattice.Crystal;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

public class CrystalHcp extends Crystal {
    
    /**
     * Cubic hcp crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public CrystalHcp() {
        this(Default.ATOM_SIZE);
    }
    
    public CrystalHcp(double latticeConstant) {
        this(new PrimitiveHexagonal(Space3D.INSTANCE));
        primitive = (PrimitiveHexagonal)((BravaisLattice)lattice).getPrimitive();
        primitive.setA(latticeConstant);
        primitive.setC(Math.sqrt(8.0/3.0)*latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private CrystalHcp(PrimitiveHexagonal primitive) {
        super(new BravaisLattice(primitive), new BasisHcp(primitive));
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