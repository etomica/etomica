/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.space.Space;
import etomica.space.Vector;
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
        setEdgeLengths(xLength, yLength);
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
    public boolean contains(etomica.space.Vector v) {
        double x = v.x(0)-position.x(0);
        double y = v.x(1)-position.x(1);
        return (x>=nX) && (x<=pX) && (y>=nY) && (y<=pY);
    }
       
    /**
     * Sets the two edge lengths of the rectangle, taking
     * each element of the given vector for the length of
     * the corresponding rectangle edge.
     */
    public void setEdgeLengths(Vector e) {
        setEdgeLengths(e.x(0), e.x(1));
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
     * @param edgeLength the new length of the sides of the square
     */
    public void setEdgeLengths(double edgeLengthX, double edgeLengthY) {
        edgeLengths.E(edgeLengthX, edgeLengthY);
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
