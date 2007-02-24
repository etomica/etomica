package etomica.math.geometry;
import etomica.exception.MethodNotImplementedException;
import etomica.space.IVector;

/**
 * Representation of a mathematical polyhedron, a 3-dimensional polytope.
 */
public class PolyhedronGeneral extends Polyhedron {
    
    protected PolyhedronGeneral(Polygon[] faces) {
        super(faces);
    }
    
    /**
     * Calculates vertices from their internal representation, which
     * for a general polyhedron does nothing because the vertices are
     * the internal representation.
     */
    public void updateVertices() {
        //does nothing, becuse in this case the vertices are the official
        //representation of the polytope
    }
    
    /**
     * Returns the 3-D volume of the polyhedron, which is its conventional volume
     */
    //must override in subclass (until a general algorithm is implemented)
    public double getVolume() {
        throw new MethodNotImplementedException();
    }
    
    /**
     * Returns true if the given point lies inside or on an edge of the polyhedron
     */
    public boolean contains(IVector vector) {
        throw new MethodNotImplementedException("General formula for 'contains' not in place");
    }

}
