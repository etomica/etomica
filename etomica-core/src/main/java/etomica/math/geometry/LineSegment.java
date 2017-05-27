/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space1d.Vector1D;
import etomica.space3d.Space3D;

/**
 * A geometrical line segment, a 1-dimensional polytope.
 * 
 * @author David Kofke
 *  
 */
public class LineSegment extends Polytope implements Rectangular {

    public LineSegment(Space embeddedSpace) {
        this(embeddedSpace, embeddedSpace.makeVector(), embeddedSpace.makeVector());
    }
    
    /**
     * Forms the segment using the given vectors as the instances used to
     * represent the end points.
     */
    LineSegment(Space embeddedSpace, Vector v0, Vector v1) {
        super(new Point[] {new Point(embeddedSpace, v0), new Point(embeddedSpace, v1)});
        edges = new LineSegment[]{this};
    }

    public void updateVertices() {
        //does nothing, becuse in this case the vertices are the representation of the polytope
    }
    
    /**
     * Returns true if the given vector lies between the ends of the segment, on
     * the line joining them. Point must lie exactly on the line to return true,
     * so rounding error will in most cases lead to return of false unless
     * embedded in a 1D space.
     */
    public boolean contains(Vector v) {
        double length = getLength();
        //(v-v0) dot (v1-v0)
        double dot = v.dot(vertices[1]) - vertices[1].dot(vertices[0])
                - v.dot(vertices[0]) + vertices[0].squared();
        if (dot < 0.0 || dot > length * length)
            return false;
        //true if (v-v0)^2 * (v1-v0)^2 equals dot^2
        return (v.Mv1Squared(vertices[0]) * (length * length) == dot * dot);
    }
    
    public double getVolume() {
        return getLength();
    }
    
    public void setVertex1(Vector newVertex1) {
        vertices[0].E(newVertex1);
    }
    
    public void setVertex2(Vector newVertex2) {
        vertices[1].E(newVertex2);
    }
    
    public double getLength() {
        return Math.sqrt(vertices[1].Mv1Squared(vertices[0]));
    }
    
    /**
     * Sets the length of the segment to the given value, keeping
     * its orientation and the position of its midpoint fixed.
     */
    public void setLength(double newLength) {
        if(vertices[0].equals(vertices[1])) {
            //currently zero length; just expand along x axis from present point
            vertices[0].setX(0,vertices[0].getX(0)-0.5*newLength);
            vertices[1].setX(0,vertices[1].getX(1)+0.5*newLength);
        } else {//not currently zero length
            double length = getLength();
            double dv = 0.5*(newLength/length - 1.0);
            //v1 <- v1 + dv*(v1-v0)
            vertices[1].TE(1.0+dv);
            vertices[1].PEa1Tv1(-dv,vertices[0]);
            //v0 <- v0 - dv*(v1-v0), but must use new v1
            vertices[0].TE(1.0+dv-dv*dv/(1.0+dv));
            vertices[0].PEa1Tv1(-dv/(1.0+dv),vertices[1]);
        }
    }

    /**
     * Sets the length equal to the element of the (presumably 1D) vector.
     * Implementation of Rectangular interface.
     */
    public void setEdgeLengths(Vector v) {
        setLength(v.getX(0));
    }
    
    /**
     * Returns the length of the line segment as the element of a 1D vector.
     * Implmentation of the Rectangular interface.
     */
    public Vector getEdgeLengths() {
        return new Vector1D(getLength());
    }
    
    public LineSegment[] getEdges() {
    	return edges;
    }


    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        LineSegment segment = new LineSegment(space, space.makeVector(new double[]{1,4}), space.makeVector(new double[]{2,8}));
        Vector p1 = space.makeVector();
        p1.setX(0, 1.5);
        p1.setX(1, 6.0);
        System.out.println(segment.contains(p1));
    }
    
    private static final long serialVersionUID = 1L;
    private final LineSegment[] edges;
}
