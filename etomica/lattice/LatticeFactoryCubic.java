package etomica.lattice;
import etomica.*;

//02.04.02

public class LatticeFactoryCubic extends BravaisLattice.Factory {
    
    /**
     * Create a square lattice (same number of sites in each dimension) of 
     * spatial dimension "D" with "dimension" elements in each direction, 
     * using a cubic primitive and occupied by the sites formed by the given factory.
     */
    public LatticeFactoryCubic(Simulation sim, AtomFactory siteFactory, 
                                int D, int dimension, double latticeConstant) {
        this(sim, siteFactory, dArray(D, dimension), latticeConstant);
    }
    /**
     * Creates a Bravais lattice with the given basis.  Number of lattice sites in each
     * dimension is determined by the given integer array.  Size of this array also
     * determines the dimension of the lattice.
     */
    public LatticeFactoryCubic(Simulation sim, AtomFactory siteFactory, 
                                int[] dimensions, double latticeConstant) {
        super(sim, siteFactory, dimensions, new PrimitiveCubic(sim));
        ((PrimitiveCubic)primitive).setSize(latticeConstant);
    }
    
    /**
     * Returns an array of length D with each element equal to dim.
     * Used to construct dimensions array to make a square lattice.
     */
    private static int[] dArray(int D, int dim) {
        int[] dimensions = new int[D];
        for(int i=0; i<D; i++) dimensions[i] = dim;
        return dimensions;
    }
}//end of LatticeFactoryCubic