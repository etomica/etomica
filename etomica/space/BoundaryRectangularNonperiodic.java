package etomica.space;

import etomica.api.IRandom;
import etomica.api.IVector;

/**
 * Boundary that is not periodic in any direction.  Volume is specified,
 * but nothing about boundary enforces the dimensions on the atoms.  This effect must
 * be introduced with a potential or other construct that is made consistent
 * with the boundary.
 */
public class BoundaryRectangularNonperiodic extends BoundaryRectangular {


    /**
     * Make a boundary with unit volume.
     */
    public BoundaryRectangularNonperiodic(Space space, IRandom random) {
        super(space, random, new boolean[space.D()], 1.0);//boolean elements will all be false
        zero = space.makeVector();
    }

    /**
     * Returns a vector with all elements zero.
     */
    public IVector centralImage(IVector r) {
        zero.E(0.0);
        return zero;
    }

    /**
     * Does nothing.
     */
    public void nearestImage(IVector dr) {
    }

    /**
     * Returns a zero-length vector.
     */
    public float[][] getOverflowShifts(IVector r, double distance) {
        return shift0;
    }

    /**
     * Returns a zero-length vector.
     */
    public double[][] imageOrigins(int nShells) {
        return origins;
    }

    private final IVector zero;
    protected final double[][] origins= new double[0][0];//cannot be static because several boxs may be using at once
    private static final long serialVersionUID = 1L;
}
