package etomica.lattice;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.space.Space;
import etomica.util.Default;

/**
 * A simple cubic lattice, with one site per cubic unit cell.
 */

public class LatticeCubicSimple extends BravaisLattice implements CubicLattice {
    
	/**
	 * Simple 3D cubic lattice with a lattice constant given by Default.ATOM_SIZE. 
	 */
    public LatticeCubicSimple() {
        this(3, Default.ATOM_SIZE);
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
        primitive.setCubicSize(latticeConstant);
	}

    /**
     * The lattice constant is the size of the cubic primitive vectors.
     */
    public void setLatticeConstant(double latticeConstant) {
        primitive.setCubicSize(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return primitive.getCubicSize();
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
    
    PrimitiveCubic primitive;//shadows superclass field
}