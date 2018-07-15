/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

//import com.fasterxml.jackson.annotation.JsonIgnore;
// import etomica.meta.annotations.IgnoreProperty;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Representation of a mathematical polytope, which is a finite region of space
 * enclosed by a finite number of hyperplanes. Subclasses include LineSegment (a
 * 1-D polytope), Polygon (a 2-D polytope), and Polyhedron (a 3-D polytope).
 * Note that the polytope may be embedded in a space of higher dimension than
 * its geometric dimension (e.g., a square -- a 2-D polytope -- by be used in a
 * three-dimensional space. A primary feature of these classes is the set of
 * points representing the vertices of the polytope. These points have the
 * dimension of the embedded space. Subclasses usually have their own internal
 * representation of the polytope that is used to calculate the vertices. It
 * restricts access to this representation to restrict the shape of the polytope
 * (e.g., a square polytope can be modified only by setting the common length of
 * its sides). A polytope can be positioned to any point in space and (to be
 * developed) given a rigid-body rotation.
 * 
 * @author kofke
 *  
 */
public abstract class Polytope implements Shape, java.io.Serializable {

    /**
     * Constructs a polytope of dimension D from the D-1 finite "planes" that
     * constitute it. Constructor does not find the vertices by calculating
     * plane intersections; instead the given planes should be configured with
     * vertices that form the polytope.
     */
    protected Polytope(Polytope[] hyperPlanes) {
        D = 1 + hyperPlanes[0].D;
        embeddedSpace = hyperPlanes[0].embeddedSpace;
        if (D > embeddedSpace.D())
            throw new IllegalArgumentException(
                    "Cannot create a polytope of dimension "
                            + D
                            + " and embed it in a lower-dimensional space (of dimension "
                            + embeddedSpace.D() + ")");
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
        this.vertices = new Vector[] { vertex };
        position = vertex;
        this.hyperPlanes = new Polytope[0];
    }

    /**
     * This method should be defined by the subclass to calculate the vertices
     * in accordance with the internal representation of the polytope. If the
     * internal representation is kept by the vertices themselves, then this
     * method need do nothing; if this method applies another internal
     * representation to the vertices (for example, it might represent a square
     * by keeping the edge length), it should then invoke
     * applyTranslationRotation() before returning. Subclass methods allowing
     * mutation of the internal representation should invoke updateVertices to
     * ensure the vertices reflect the internal state. Vertices may be
     * changed directly, but these changes will overridden by updateVertices if
     * it is called afterwards.
     */
    public abstract void updateVertices();

    /**
     * Calculate transformed vertices from current values of vertices, position,
     * and orientation. Does not invoke updateVertices before transforming, and
     * so should be called only by subclass and only at the end of the updateVertices
     * method. Otherwise translation and rotation will be performed from current
     * value of vertices, compounding previous application of
     * translation/rotation. NOTE: Rotation is not yet implemented.
     */
    protected void applyTranslationRotation() {
        if(noTranslation) return;
        for (int i = 0; i < vertices.length; i++) {
            vertices[i].Ev1Pv2(position, vertices[i]);
        }
    }

    /**
     * Returns an array of all vertices of the polytrope. For most subclasses
     * this array will not provide the direct representation of the polytope.
     * Instead polytope will be defined in another internal representation, which is
     * used to calculate the vertices via the updateVertices method. In a few
     * subclasses the vertices <i>are</i> the representation, so alteration of the
     * returned array will change the polytope.
     */
    public Vector[] getVertices() {
        return vertices;
    }

    /**
     * Sets the position of the geometric center of the polytope.
     */
    public void setPosition(Vector r) {
        position.E(r);
        noTranslation = position.isZero();
        updateVertices();
    }

    /**
     * Returns the geometric center of the polytope. The returned vector is used
     * to define the position internally, so changes to it will affect the
     * logical position of the polytope, without immediately changing the vertices.
     * Use setPosition to change the polytope position.
     */
    public Vector getPosition() {
        return position;
    }

    /**
     * Returns the (hyper)volume enclosed by the polytope.  So
     * for a LineSegment instance, this is its length; for a
     * polygon, it is the area; for a polyhedron, it is the volume.  
     */
    public abstract double getVolume();

    /**
     * Returns <code>true</code> if the given vector lies inside the polytope,
     * <code>false</code> otherwise.
     */
    public abstract boolean contains(Vector v);

    /**
     * Number of vertices in the polytrope. A vertex is a point where D edges
     * meet.
     */
//    @JsonIgnore
    public final int getVertexCount() {
        return vertices.length;
    }

    /**
     * The space defining the vectors used to represent the vertices.
     */
//    @JsonIgnore
    public Space getEmbeddedSpace() {
        return embeddedSpace;
    }

    /**
     * The dimension of this polytope, which is not necessarily the dimension
     * of the space it is embedded in.  For example, for a cube D = 3, and for
     * a square D = 2.
     */
    public int D() {
        return D;
    }

    /**
     * Returns a list of all vertices appearing in the given array of
     * polytopes.  Each vertex appears in the list only once.  Used
     * by constructor.
     */
    private static Vector[] allVertices(Polytope[] hyperPlanes) {
        List<Vector> list = new ArrayList<>();
        for (int i = 0; i < hyperPlanes.length; i++) {
            Vector[] vertices = hyperPlanes[i].vertices;
            for (int j = 0; j < vertices.length; j++) {
                if (!list.contains(vertices[j])) {
                    list.add(vertices[j]);
                }
            }
        }
        return list.toArray(new Vector[0]);
    }
    
    public abstract LineSegment[] getEdges();
    
    public String toString() {
        StringBuilder str = new StringBuilder();
        str.append("{");
        for(int i=0; i<hyperPlanes.length; i++) {
            str.append(hyperPlanes[i].toString());
            if(i < hyperPlanes.length-1) str.append(",");
        }
        str.append("}");
        return str.toString();
    }

    protected final Space embeddedSpace;
    protected final int D;

    /**
     * These vertices are used in all external representations of the polygon.
     * Changes to them do not necessarily change the state of the polygon.
     */
    protected final Vector[] vertices;
    protected final Vector position;
    protected final Polytope[] hyperPlanes;
    private boolean noTranslation = true;
    //    protected final Orientation orientation;

}
