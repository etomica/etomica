package etomica.lattice;
import etomica.Space;

public class LatticeCubic extends Lattice {
    
    /**
     * Create a cubic Bravais lattice (same number of sites in each dimension) of dimension "d" 
     * with "size" elements in each dimension, occupied by the given basis
     */
    public LatticeCubic(int d, int size, double pVectorLength, Basis basis) {
        super(d, size, pVectors(pVectorLength, d), basis);
    }
    /**
     * Creates a Bravais lattice with the given basis.  Number of lattice sites in each
     * dimension is determined by the given integer array.  Size of this array also
     * determines the dimension of the lattice.
     */
    public LatticeCubic(int[] dimensions, double pVectorLength, Basis basis) {
        super(dimensions, pVectors(pVectorLength, dimensions.length), basis);
    }
    
    /**
     * Makes an orthogonal set of vectors, each of the given length, for a d-dimensional space.
     */
    private static Space.Vector[] pVectors(double length, int d) {
        Space.Vector[] vectors = new Space.Vector[d];
        for(int i=0; i<d; i++) {
            vectors[i] = Space.makeVector(d);
            vectors[i].setComponent(i,length);
        }
        return vectors;
    }
}