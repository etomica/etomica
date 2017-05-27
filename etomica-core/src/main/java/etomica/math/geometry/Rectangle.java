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
 * A 2-dimensional geometric rectangle.
 * @author kofke
 *
 */
public class Rectangle extends Polygon implements Rectangular {

    /**
     * Constructs a square rectangle of unit size.
     */
    public Rectangle(Space embeddedSpace) {
        this(embeddedSpace, 1.0, 1.0);
    }
    
    /**
     * Constructs a rectangle with edge lengths having the given values.
     */
    public Rectangle(Space embeddedSpace, double xLength, double yLength) {
        super(embeddedSpace, 4);
        setEdgeLengths(embeddedSpace.makeVector(new double[]{xLength, yLength}));
    }

    /**
     * Returns the area enclosed by the rectangle.
     */
    public double getArea() {
        return 4 * pX * pY;
    }
    
    /**
     * Returns the perimeter of the rectangle.
     */
    public double getPerimeter() {
        return 2 * 2 * (pX + pY);
    }

    
    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean contains(Vector v) {
        double x = v.getX(0)-position.getX(0);
        double y = v.getX(1)-position.getX(1);
        return (x>=nX) && (x<=pX) && (y>=nY) && (y<=pY);
    }
       
    /**
     * Sets the two edge lengths of the rectangle, taking
     * each element of the given vector for the length of
     * the corresponding rectangle edge.
     */
    public void setEdgeLengths(Vector e) {
        edgeLengths.E(e);
        setEdgeLengths(e.getX(0), e.getX(1));
    }
    
    /**
     * Returns a vector with elements equal to the
     * edge lengths of the rectangle.  The returned vector is not
     * used to represent the rectangle internally, so changing its
     * values will not affect the state of the rectangle.
     */
    public Vector getEdgeLengths() {
        return edgeLengths;
    }

    /**
     * @param edgeLengthX the new length of the sides of the square
     * @param edgeLengthY the new length of the sides of the square
     */
    protected void setEdgeLengths(double edgeLengthX, double edgeLengthY) {
        nX = -0.5*edgeLengthX;
        pX = +0.5*edgeLengthX;
        nY = -0.5*edgeLengthY;
        pY = +0.5*edgeLengthY;
        updateVertices();
    }
    
    public void updateVertices() {
        ((Vector2D)vertices[0]).E(nX,nY);
        ((Vector2D)vertices[1]).E(nX,pY);
        ((Vector2D)vertices[2]).E(pX,pY);
        ((Vector2D)vertices[3]).E(pX,nY);
        applyTranslationRotation();
    }
    
    private final Vector2D edgeLengths = new Vector2D();
    private double nX, pX;//n = -size/2, p = +size/2
    private double nY, pY;//n = -size/2, p = +size/2


}
