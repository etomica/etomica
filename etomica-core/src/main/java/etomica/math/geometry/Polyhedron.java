/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

//import com.fasterxml.jackson.annotation.JsonIgnore;
import etomica.space.Vector;

/**
 * Representation of a mathematical polyhedron, a 3-dimensional polytope. Contains
 * all information needed to represent a polyhedron, methods to set its position
 * and orientation, and methods that return an external representation of it. 
 * Provides no external means to change the polyhedron size and shape, as the form
 * and capabilities of mutator methods are particular to the type of polyhedron
 * defined by the subclass. 
 */
public abstract class Polyhedron extends Polytope {
    
    /**
     * Constructs a polyhedron with the given faces for its sides.  Faces
     * should be constructed so they have the correct sharing of edges
     * between adjacent faces, and sharing of vertices between intersecting
     * edges.
     */
    protected Polyhedron(Polygon[] faces) {
        super(faces);
        this.faces = faces;
        this.edges = allEdges(faces);
    }
    
    /**
     * Returns all faces defined by the polyhedron.
     */
//    @JsonIgnore
    public Polygon[] getFaces() {
        return faces;
    }
    
    /**
     * Returns all edges defined by the polyhedron.
     */
    public LineSegment[] getEdges() {
        return edges;
    }
        
    /**
     * Returns the sum of the length of the edges of the polyhedron
     */
    public double getPerimeter() {
        double sum = 0.0;
        for(int i=0; i<edges.length; i++) {
            sum += edges[i].getLength();
        }
        return sum;
    }

    /**
     * Returns the sum the area of the faces of the polyhedron
     */
    public double getSurfaceArea() {
        double sum = 0.0;
        for(int i=0; i<faces.length; i++) {
            sum += faces[i].getArea();
        }
        return sum;
    }
    
    /**
     * Returns the perpendicular distance to the nearest face of the 
     * polyhedron.  Assumes that given point is contained in polyhedron;
     * if not, returns NaN.
     * @param r
     * @return
     */
    public double distanceTo(Vector r) {
        if(!contains(r)) return Double.NaN;
        double d = Double.POSITIVE_INFINITY;
        Plane plane = new Plane(embeddedSpace);
        for(int i=0; i<faces.length; i++) {
            Polygon f = faces[i];
            plane.setThreePoints(vertices[0], 
                    vertices[1], vertices[2]);
            double d1 = Math.abs(plane.distanceTo(r));
            if(d1 < d) d = d1;
        }
        return d;
    }
    
    /**
     * Finds all edge instances in the given array of faces, and returns
     * them in an array.  Each edge appears in the array once, although it
     * should be part of multiple faces.  Used by constructor.
     */
    private static LineSegment[] allEdges(Polygon[] faces) {
        List<LineSegment> list = new ArrayList<>();
        for(int i=0; i<faces.length; i++) {
            LineSegment[] edges= faces[i].getEdges();
            for(int j=0; j<edges.length; j++) {
                if(!list.contains(edges[j])) {
                    list.add(edges[j]);
                }
            }
        }
        return list.toArray(new LineSegment[0]);
    }

    protected final LineSegment[] edges;
    protected final Polygon[] faces;
    
 }
