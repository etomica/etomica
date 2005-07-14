package etomica.space;

import etomica.Space;


/**
 * Boundary that is not periodic in any direction.  Volume is specified,
 * but nothing about boundary enforces dimensions on the atoms.  This must
 * be introduced with a potential or other construct that is made consistent
 * with the boundary.
 */

/*
 * History
 * Created on Jan 27, 2005 by kofke
 */
public class BoundaryRectangularNonperiodic extends BoundaryRectangular {

    /**
     * @param space
     */
    public BoundaryRectangularNonperiodic(Space space) {
        super(space, new boolean[space.D()]);//boolean elements will all be false
        zero = space.makeVector();
    }

    /**
     * Returns a vector with all elements zero.
     */
    public Vector centralImage(Vector r) {
        zero.E(0.0);
        return zero;
    }

    /**
     * Does nothing.
     */
    public void nearestImage(Vector dr) {
    }

    /**
     * Returns a zero-length vector.
     */
    public float[][] getOverflowShifts(Vector r, double distance) {
        return shift0;
    }

    /**
     * Returns a zero-length vector.
     */
    public double[][] imageOrigins(int nShells) {
        return origins;
    }

    private final Vector zero;
    protected final float[][] shift0 = new float[0][0];//cannot be static because several phases may be using at once
    protected final double[][] origins= new double[0][0];//cannot be static because several phases may be using at once
    
}
