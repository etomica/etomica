/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IBoundary;
import etomica.api.IBoundaryEvent;
import etomica.api.IBoundaryEventManager;
import etomica.api.IBox;
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
public abstract class Boundary implements IBoundary, java.io.Serializable {

    /**
     * Subclasses must invoke this constructor and provide a Space instance that
     * can be used to generate Vectors, and a Polytope that defines the shape
     * and volume of the boundary region. Both are final.
     */
    public Boundary(ISpace space, Polytope shape) {
        this.space = space;
        this.shape = shape;
        double zip[] = new double[space.D()];
        for(int i = 0; i < space.D(); i++) zip[i] = 0.0;
        center = space.makeVector(zip);
        eventManager = new BoundaryEventManager();
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
            inflateEvent = new BoundaryEvent(this);
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

    public IVector getCenter() {
        return center;
    }
    
    public IBoundaryEventManager getEventManager() {
        return eventManager;
    }
    
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
    private final IVector center;
    protected final Polytope shape;
    protected final ISpace space;
    protected IBox box;
    protected IBoundaryEvent inflateEvent;
    protected BoundaryEventManager eventManager;
}
