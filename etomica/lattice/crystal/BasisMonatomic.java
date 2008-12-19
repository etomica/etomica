package etomica.lattice.crystal;

import etomica.api.IVectorMutable;
import etomica.space.ISpace;

/**
 * Single-atom basis with the coordinate at the origin.
 */
public class BasisMonatomic extends Basis {
    
    /**
     * Creates a single-atom basis with the coordinate at the origin.
     */
    public BasisMonatomic(ISpace space) {
        super(new IVectorMutable[] {space.makeVector()});
    }
    	
    private static final long serialVersionUID = 1L;
}
