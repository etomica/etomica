/*
 * History
 * Created on Nov 25, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.Space;
import etomica.Space.Vector;

/**
 * Wraps a Polytope so that it may be given a position and 
 * orientation different from the fixed configuration defined
 * for the elementary polytopes.
 * @author kofke
 *
 */
public class PolytopeMovable extends Polytope {

    /**
     * @param space
     */
    public PolytopeMovable(Polytope polytope) {
        super(polytope.space());
        this.polytope = polytope;
        position = space.makeVector();
        orientation = space.makeOrientation();
        throw new etomica.exception.MethodNotImplementedException("Class still under construction");
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Polytope#setSize(double)
     */
    public void setSize(double size) {
        polytope.setSize(size);
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Polytope#getSize()
     */
    public double getSize() {
        return polytope.getSize();
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Polytope#volume()
     */
    public double volume() {
        return polytope.volume();
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Polytope#vertex()
     */
    public Vector[] vertex() {
        //get vertices from wrapped polytope and translate/rotate them
        return null;
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Polytope#inCell(etomica.Space.Vector)
     */
    public boolean inCell(Vector v) {
        //do rotations/translations and determine inCell
        return false;
    }
    
    /**
     * @return the wrapped polytope.
     */
    public Polytope getPolytope() {
        return polytope;
    }
    
    private final Polytope polytope;
    private final Space.Vector position;
    private final Space.Orientation orientation;
}
