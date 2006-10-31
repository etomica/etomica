package etomica.lattice.crystal;

import etomica.space.Space;
import etomica.space.Vector;

/**
 * Single-atom basis with the coordinate at the origin.
 */
public class BasisMonatomic extends Basis {
    
    /**
     * Creates a single-atom basis with the coordinate at the origin.
     */
    public BasisMonatomic(Space space) {
        super(new Vector[] {space.makeVector()});
    }
    	
    private static final long serialVersionUID = 1L;
}
