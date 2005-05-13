/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import java.util.LinkedList;

import etomica.Space;
import etomica.space.Vector;

/**
 * Representation of a mathematical polytope, which is a finite region
 * of space enclosed by a finite number of hyperplanes.  Subclasses
 * include LineSegment (a 1-D polytope), Polygon (a 2-D polytope), and
 * Polyhedron (a 3-D polytope).
 * 
 * @author kofke
 *
 */
public abstract class Polytope {

    /**
     * Constructs a polytope of dimension D from the D-1 "planes"
     * that constitute it.
     */
    protected Polytope(Polytope[] hyperPlanes) {
        D = 1 + hyperPlanes[0].D;
        embeddedSpace = hyperPlanes[0].embeddedSpace;
        if(D > embeddedSpace.D()) throw new IllegalArgumentException("Cannot create a polytope of dimension "+D+" and embed it in a lower-dimensional space (of dimension "+embeddedSpace.D()+")");
        this.vertices = allVertices(hyperPlanes);
        position = embeddedSpace.makeVector();
        this.hyperPlanes = hyperPlanes;
    }
    
    /**
     * Constructor used for the Point subclass
     */
    Polytope(Space embeddedSpace, Vector vertex) {
        D = 0;
        this.embeddedSpace = embeddedSpace;
        this.vertices = new Vector[] {vertex};
        position = vertex;
        this.hyperPlanes = new Polytope[0];
    }
    
    /**
     * This method should be defined by the subclass to calculate
     * the vertices in accordance with the internal representation
     * of the polytope.  If the internal representation is kept by
     * the vertices themselves, then this method need do nothing;
     * if this method applies another internal representation to the
     * vertices (for example, it might represent a square by keeping
     * the edge length), it should then invoke 
     * applyTranslationRotation() before returning.
     */
    public abstract void updateVertices();
    
    /**
     * Calculated transformed vertices from current values of vertices,
     * position, and orientation. Modification of the untransformed vertices
     * is performed by methods defined by the subclass.
     */
    protected void applyTranslationRotation() {
        for (int i = 0; i < vertices.length; i++) {
            vertices[i].Ev1Pv2(position, vertices[i]);
        }
    }

    /**
     * Returns an array of all vertices of the polytrope.
     */
    public Vector[] vertices() {
        updateVertices();
        return vertices;
    }

     /**
     * Sets the position of the geometric center of the polytope.
     */
    public void setPosition(Vector r) {
        position.E(r);
        updateVertices();
    }

    /**
     * Returns the geometric center of the polytope. The returned vector is used
     * to define the position internally, so changes to it will affect the
     * position of the polytope.
     */
    public Vector getPosition() {
        return position;
    }
    
    public abstract double getVolume();
    
    /**
     * Returns <code>true</code> if the given vector lies inside the 
     * polytope, <code>false</code> otherwise.
     */
    public abstract boolean contains(Vector v);

    /**
     * Number of vertices in the polytrope.
     * A vertex is a point where D edges meet.
     */
    public final int vertexCount() {return vertices.length;}
    
    public Space embeddedSpace() {
        return embeddedSpace;
    }
    
    public int D() {
        return D;
    }

    private static Vector[] allVertices(Polytope[] hyperPlanes) {
        LinkedList list = new LinkedList();
        for(int i=0; i<hyperPlanes.length; i++) {
            Vector[] vertices = hyperPlanes[i].vertices;
            for(int j=0; j<vertices.length; j++) {
                if(!list.contains(vertices[j])) {
                    list.add(vertices[j]);
                }
            }
        }
        return (Vector[])list.toArray(new Vector[0]);
    }
    
    protected final Space embeddedSpace;
    protected final int D;
    
    /**
     * These vertices are used internally to represent the state (size, shape)
     * of the polygon. They are defined with respect to the center of the
     * polygon as it lies in the x-y plane. The subclass is responsible for
     * defining and updating these vertices.
     */
 //   protected final Vector[] vertices;
    /**
     * These vertices are used in all external representations of the polygon.
     * Changes to them do not change the state of the polygon.
     */
    protected final Vector[] vertices;
    protected final Vector position;
    protected final Polytope[] hyperPlanes;
//    protected final Orientation orientation;

}
