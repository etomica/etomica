package etomica.space;

import etomica.Default;
import etomica.Space;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 27, 2005 by kofke
 */
public class BoundaryNone extends Boundary {

    private final Vector dimensions;
    
    /**
     * @param space
     */
    public BoundaryNone(Space space) {
        super(space);
        dimensions = space.makeVector();
        zero = space.makeVector();
        temp = space.makeVector();
        dimensions.E(Default.BOX_SIZE);
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#centralImage(etomica.space.Vector)
     */
    public Vector centralImage(Vector r) {
        zero.E(0.0);
        return zero;
    }

    /* (non-Javadoc)
     * @see etomica.NearestImageTransformer#nearestImage(etomica.space.Vector)
     */
    public void nearestImage(Vector dr) {
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#volume()
     */
    public double volume() {
        return dimensions.productOfElements();
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#dimensions()
     */
    public Vector dimensions() {
        return dimensions;
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#setDimensions(etomica.space.Vector)
     */
    public void setDimensions(Vector v) {
        dimensions.E(v);
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#randomPosition()
     */
    public Vector randomPosition() {
        temp.setRandomCube();
        temp.PE(0.5);
        temp.TE(dimensions);
        return temp;
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#getOverflowShifts(etomica.space.Vector, double)
     */
    public float[][] getOverflowShifts(Vector r, double distance) {
        return shift0;
    }

    /* (non-Javadoc)
     * @see etomica.space.Boundary#imageOrigins(int)
     */
    public double[][] imageOrigins(int nShells) {
        return origins;
    }

    private final Vector zero;
    private final Vector temp;
    protected final float[][] shift0 = new float[0][0];//cannot be static because several phases may be using at once
    protected final double[][] origins= new double[0][0];//cannot be static because several phases may be using at once
    
}
