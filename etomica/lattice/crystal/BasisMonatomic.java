package etomica.lattice.crystal;

import etomica.Space;
import etomica.lattice.Basis;
import etomica.space.Vector;

/**
 * Single-atom basis with the coordinate at the origin.
 */
 
public class BasisMonatomic implements Basis, java.io.Serializable {
    
    /**
     * Creates a single-atom basis with the coordinate at the origin.
     * @param D the spatial dimension of the crystal
     */
    public BasisMonatomic(Space space) {
        coordinates = new Vector[] {space.makeVector()};
    }
    	
	/**
	 * Number of atoms in the basis.
	 * @return int
	 */
	public int size() {
        return 1;
    }
    
    /**
     * Returns an array with a single vector at the origin.
     */
    public Vector[] positions() {
        return coordinates;
    }

    private Vector[] coordinates;
}
