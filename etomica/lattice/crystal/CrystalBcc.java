package etomica.lattice.crystal;
import etomica.Default;
import etomica.Space3D;
import etomica.lattice.BravaisLattice;
import etomica.lattice.Crystal;

/**
 * Cubic primitive with a 2-site bcc basis.
 */

public class CrystalBcc extends Crystal {
    
	/**
	 * Cubic bcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of size Default.ATOM_SIZE. 
	 */
    public CrystalBcc() {
        this(2.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
    }
    
	public CrystalBcc(double latticeConstant) {
		this(new PrimitiveCubic(Space3D.INSTANCE));
        primitive = (PrimitiveCubic)((BravaisLattice)lattice).getPrimitive();
        primitive.setSize(latticeConstant);
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
	 * new BasisCubicBcc (which needs the new primitive) to super.
	 */	
	private CrystalBcc(PrimitiveCubic primitive) {
		super(new BravaisLattice(primitive), new BasisCubicBcc(primitive));
	}
    
    /**
     * Returns the primitive the determines the lattice constant.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveCubic primitive() {
        return primitive;
    }

    /**
     * Returns "Bcc".
     */
    public String toString() {return "Bcc";}
    
    private PrimitiveCubic primitive;
    
}//end of CrystalBcc