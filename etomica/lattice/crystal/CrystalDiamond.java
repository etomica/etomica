package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.BravaisLattice;
import etomica.lattice.Crystal;

/**
 * Cubic primitive with a 4-site fcc basis, on which each site 
 * is a 2-site diamond basis.
 */

 /* History
  * 09/26/02 (DAK) new
  * 01/20/04 (DAK) revised constructors; added one taking atomFactory argument
  */
public class CrystalDiamond extends Crystal {
    
    /**
     * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public CrystalDiamond() {
        this(4.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
    }
    
    public CrystalDiamond(double latticeConstant) {
        this(new PrimitiveCubic(Space3D.INSTANCE));
        primitive = (PrimitiveCubic)((BravaisLattice)lattice).getPrimitive();
        primitive.setSize(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private CrystalDiamond(PrimitiveCubic primitive) {
        super(new BravaisLattice(primitive), new BasisCubicDiamond(primitive));
    }
    
    /**
     * Returns the primitive the determines the lattice constant.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveCubic primitive() {
        return primitive;
    }

    /**
     * Returns "Diamond".
     */
    public String toString() {return "Diamond";}
    
    private PrimitiveCubic primitive;
    
}//end of CrystalDiamond