package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.BravaisLattice;
import etomica.lattice.Crystal;

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
public class CrystalFcc extends Crystal {

    /**
     * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
     */
    public CrystalFcc() {
        this(Math.sqrt(2.0)*Default.ATOM_SIZE);
    }
    
    public CrystalFcc(double latticeConstant) {
        this(new PrimitiveCubic(Space3D.INSTANCE));
        primitive = (PrimitiveCubic)((BravaisLattice)lattice).getPrimitive();
        primitive.setSize(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicBcc (which needs the new primitive) to super.
     */ 
    private CrystalFcc(PrimitiveCubic primitive) {
        super(new BravaisLattice(primitive), new BasisCubicFcc(primitive));
    }
    
    /**
     * Returns the primitive the determines the lattice constant.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveCubic primitive() {
        return primitive;
    }

    /**
     * Returns "Fcc".
     */
    public String toString() {return "Fcc";}
    
    private PrimitiveCubic primitive;
}//end of CrystalFcc