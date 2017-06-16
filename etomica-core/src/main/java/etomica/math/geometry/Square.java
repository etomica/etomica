/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space2d.Vector2D;

/**
 * A 2-dimensional geometric square.
 * @author kofke
 *
 */
public class Square extends Polygon {

    /**
     * Constructs a square of unit size.
     */
    public Square(Space embeddedSpace) {
        this(embeddedSpace, 1.0);
    }
    
    /**
     * Constructs a square with edge length having the given value.
     * @param size edge length of the cube
     */
    public Square(Space embeddedSpace, double size) {
        super(embeddedSpace, 4);
        setEdgeLength(size);
    }

    /**
     * Returns edgeLength^2.
     */
    public double getArea() {
        return edgeLength*edgeLength;
    }

    
    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean contains(Vector v) {
        double x = v.getX(0)-position.getX(0);
        double y = v.getX(1)-position.getX(1);
        return (x>=n) && (x<=p) && (y>=n) && (y<=p);
    }
       
    /**
     * @return the length the edge of the square.
     */
    public double getEdgeLength() {
        return edgeLength;
    }
    
    /**
     * @param edgeLength the new length of the sides of the square
     */
    public void setEdgeLength(double edgeLength) {
        this.edgeLength = edgeLength;
        n = -0.5*edgeLength;
        p = +0.5*edgeLength;
        updateVertices();
    }
    
    public void updateVertices() {
        ((Vector2D)vertices[0]).E(n,n);
        ((Vector2D)vertices[1]).E(n,p);
        ((Vector2D)vertices[2]).E(p,p);
        ((Vector2D)vertices[3]).E(p,n);
        applyTranslationRotation();
    }
    
    private double edgeLength;
    private double n, p;//n = -size/2, p = +size/2


}
