package etomica.lattice;
import etomica.Default;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.space3d.Space3D;

/**
 * Cubic primitive with a 4-site fcc basis.
 * 
 */

 /* History
  * 09/23/02 (DAK) new
  * 12/09/02 (skkwak) built another contructor to pass an instance of
  * AtomFactory class into BasisCubicFcc to change atoms to sites (concrete
  * sites under basis from BravisLattice)
  * 01/20/04 (DAK) restructured constructors
  */
public class LatticeCubicFcc extends LatticeCrystal implements CubicLattice {

    /**
     * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public LatticeCubicFcc() {
        this(Math.sqrt(2.0)*Default.ATOM_SIZE);
    }
    
    public LatticeCubicFcc(double latticeConstant) {
        this(new PrimitiveCubic(Space3D.INSTANCE));
        primitive = (PrimitiveCubic)crystal.getLattice().getPrimitive();
        primitive.setCubicSize(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private LatticeCubicFcc(PrimitiveCubic primitive) {
        super(new Crystal(primitive, new BasisCubicFcc(primitive)));
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
        primitive.setCubicSize(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return primitive.getCubicSize();
    }


    /**
     * Returns "Fcc".
     */
    public String toString() {return "Fcc";}
    
    private PrimitiveCubic primitive;
}//end of CrystalFcc