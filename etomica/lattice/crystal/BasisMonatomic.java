package etomica.lattice.crystal;

import etomica.api.IVector;
import etomica.space.Space;

/**
 * Single-atom basis with the coordinate at the origin.
 */
public class BasisMonatomic extends Basis {
    
    /**
     * Creates a single-atom basis with the coordinate at the origin.
     */
    public BasisMonatomic(Space space) {
        super(new IVector[] {space.makeVector()});
    }
    	
    private static final long serialVersionUID = 1L;
}
