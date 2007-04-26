package etomica.lattice;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.space.Space;

/**
 * A simple cubic lattice, with one site per cubic unit cell.
 */

public class LatticeCubicSimple extends BravaisLattice implements CubicLattice {
    
	/**
	 * Simple 3D cubic lattice with a unit lattice constant. 
	 */
    public LatticeCubicSimple() {
        this(3, 1.0);
    }
    
    /**
     * Makes a simple cubic lattice of a spatial dimension given by the first
     * argument, and with the given lattice constant.
     * @param D dimension of space of lattice; values of 2 or 3 are expected, D = 1 might also work
     * @param latticeConstant spacing between adjacent lattice sites
     */
	public LatticeCubicSimple(int D, double latticeConstant) {
		super(new PrimitiveCubic(Space.getInstance(D)));
        this.primitive = (PrimitiveCubic)getPrimitive();
        primitive.setSizeABC(latticeConstant);
	}

    /**
     * The lattice constant is the size of the cubic primitive vectors.
     */
    public void setLatticeConstant(double latticeConstant) {
        primitive.setSizeABC(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return primitive.getSizeABC();
    }
    
    /**
     * Rescales the lattice by the given factor. Multiplies the lattice constant
     * by the given value.
     */
    public void scaleBy(double scaleFactor) {
        primitive.scaleSize(scaleFactor);
    }
    
    /**
     * Returns 1.
     */
    public int getBasisSize() {
        return 1;
    }

    /**
     * Returns "Simple cubic".
     */
    public String toString() {return "Simple cubic";}
    
    private PrimitiveCubic primitive;//shadows superclass field
}