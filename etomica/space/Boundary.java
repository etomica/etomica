package etomica.space;

import etomica.api.IBoundary;
import etomica.api.IVector;
import etomica.lattice.IndexIteratorSizable;
import etomica.math.geometry.Polytope;

/**
 * Parent class of boundary objects that describe the size and periodic nature
 * of the box boundaries. Each Box has its own instance of this class. It
 * may be referenced by a coordinate pair when computing distances between
 * atoms, or by a cell iterator when generating neighbor lists. It is also used
 * by objects that enforce periodic images.
 * 
 */
public abstract class Boundary implements NearestImageTransformer, java.io.Serializable, IBoundary {

    /**
     * Subclasses must invoke this constructor and provide a Space instance that
     * can be used to generate Vectors, and a Polytope that defines the shape
     * and volume of the boundary region. Both are final.
     */
    public Boundary(Space space, Polytope shape) {
        this.space = space;
        this.shape = shape;
    }

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getShape()
	 */
    public Polytope getShape() {
        return shape;
    }

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#volume()
	 */
    public double volume() {
        return shape.getVolume();
    }

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#centralImage(etomica.api.IVector)
	 */
    public abstract IVector centralImage(IVector r);

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#nearestImage(etomica.api.IVector)
	 */
    public abstract void nearestImage(IVector dr);

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getDimensions()
	 */
    public abstract IVector getDimensions();

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#setDimensions(etomica.api.IVector)
	 */
    public abstract void setDimensions(IVector v);

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#randomPosition()
	 */
    public abstract IVector randomPosition();

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getPeriodicVectors()
	 */
    public abstract IVector[] getPeriodicVectors();
    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getIndexIterator()
	 */
    public abstract IndexIteratorSizable getIndexIterator();
    
    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getOverflowShifts(etomica.api.IVector, double)
	 */
    public abstract float[][] getOverflowShifts(IVector r, double distance);

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#imageOrigins(int)
	 */
    public abstract double[][] imageOrigins(int nShells);

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getBoundingBox()
	 */
    public abstract IVector getBoundingBox();
    
//    protected final Space space;
    protected final Polytope shape;
    protected final Space space;

}
