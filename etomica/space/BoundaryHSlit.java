package etomica.space;

import etomica.Space;

/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
/**
 * Class for implementing slit periodic boundary conditions 
 */
public class BoundaryHSlit extends BoundaryPeriodicSquare {
    
    /**
     * Constructor for periodic boundary conditions with a slit 
     * in the given dimension.
     * @param space
     * @param slitDim slit dimension (in which PBC are not imposed).
     */
    public BoundaryHSlit(Space space, int slitDim) {
        super(space);
        sDim = slitDim;
     }
    public void nearestImage(Vector dr) {
        double x = dr.x(sDim);
        super.nearestImage(dr);
        dr.setX(sDim,x);
    }
    public Vector centralImage(Vector r) {
        //superclass sets modShift
        super.centralImage(r);
        modShift.setX(sDim,0.0);
        return modShift;
    }

    private final int sDim;
}
