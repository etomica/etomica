package etomica.lattice;
import etomica.Default;
import etomica.Space3D;
import etomica.lattice.crystal.BasisCubicDiamond;
import etomica.lattice.crystal.PrimitiveCubic;

/**
 * Cubic primitive with a 4-site fcc basis, on which each site 
 * is a 2-site diamond basis.
 */

 /* History
  * 09/26/02 (DAK) new
  * 01/20/04 (DAK) revised constructors; added one taking atomFactory argument
  */
public class LatticeCubicDiamond extends LatticeCrystal implements LatticeCubic {
    
    /**
     * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public LatticeCubicDiamond() {
        this(4.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
    }
    
    public LatticeCubicDiamond(double latticeConstant) {
        this(new PrimitiveCubic(Space3D.INSTANCE));
        primitive = (PrimitiveCubic)((BravaisLattice)crystal.getLattice()).getPrimitive();
        primitive.setSize(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private LatticeCubicDiamond(PrimitiveCubic primitive) {
        super(new Crystal(new BravaisLattice(primitive), new BasisCubicDiamond(primitive)));
    }
    
    /**
     * Returns the primitive the determines the lattice constant.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveCubic primitive() {
        return primitive;
    }
    
    /**
     * The lattice constant is the size of the cubic primitive vectors
     * of the lattice underlying this crystal.
     */
    public void setLatticeConstant(double latticeConstant) {
        primitive.setSize(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return primitive.getSize();
    }
    

    /**
     * Returns "Diamond".
     */
    public String toString() {return "Diamond";}
    
    private PrimitiveCubic primitive;
    
}//end of CrystalDiamond