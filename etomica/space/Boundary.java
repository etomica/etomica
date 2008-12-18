package etomica.space;

import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IBoxEvent;
import etomica.api.IVector;
import etomica.box.BoxInflateEvent;
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
public abstract class Boundary implements IBoundary, java.io.Serializable {

    /**
     * Subclasses must invoke this constructor and provide a Space instance that
     * can be used to generate Vectors, and a Polytope that defines the shape
     * and volume of the boundary region. Both are final.
     */
    public Boundary(ISpace space, Polytope shape) {
        this.space = space;
        this.shape = shape;
    }

    /**
     * Sets the boundary's box to be the given box.  The box may be null to
     * indicate that the boundary is not associated with a box.
     */
    public void setBox(IBox newBox) {
        box = newBox;
        if (box == null) {
            inflateEvent = null;
        }
        else {
            inflateEvent = new BoxInflateEvent(box);
        }
    }

    /**
     * Returns the boundary's box.  May be null if the boundary is not
     * associated with a box.
     */
    public IBox getBox() {
        return box;
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
	 * @see etomica.space.IBoundary#randomPosition()
	 */
    public abstract IVector randomPosition();

    /* (non-Javadoc)
	 * @see etomica.space.IBoundary#getIndexIterator()
	 */
    public abstract IndexIteratorSizable getIndexIterator();
    
	/**
	 * Set of vectors describing the displacements needed to translate the
	 * central image to all of the periodic images. The first index specifies
	 * each periodic image, while the second index indicates the xyz components
	 * of the translation vector.
	 * 
	 * @param nShells
	 *            the number of shells of images to be computed
	 */
    public abstract double[][] imageOrigins(int nShells);

    private static final long serialVersionUID = 1L;
//    protected final Space space;
    protected final Polytope shape;
    protected final ISpace space;
    protected IBox box;
    protected IBoxEvent inflateEvent;
}
