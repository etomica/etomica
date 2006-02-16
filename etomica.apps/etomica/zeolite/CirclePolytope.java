package testing;

import etomica.space.Space;
import etomica.space2d.Vector2D;
import etomica.math.geometry.*;

/**
 * A 2-dimensional geometric circle that is consistent with a Polygon.
 * @author kofke
 *
 */

public class CirclePolytope extends Polygon{


    /**
     * Constructs a square of unit size.
     */
    public CirclePolytope(Space embeddedSpace) {
        this(embeddedSpace, 1.0,10);
    }
    
    /**
     * Constructs a square with edge length having the given radius.
     * @param radius of the circle
     */
    public CirclePolytope(Space embeddedSpace, double size,int sides) {
        super(embeddedSpace, sides);
        setEdgeLength(2*Math.PI*size/sides);
        radius = size;
        nsides = sides;
    }

    /**
     * Returns edgeLength^2.
     */
    public double getArea() {
        return radius*radius*Math.PI;
    }

    
    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean contains(etomica.space.Vector v) {
        double x = v.x(0)-position.x(0);
        double y = v.x(1)-position.x(1);
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
        updateVertices();
    }
    
    public void updateVertices() {
    	double degreeIncr = 360/nsides;
    	double degree = 0;
    	for(int i=0;i<nsides;i++){
    		x = radius*Math.cos(degree);
    		y = radius*Math.sin(degree);
    		degree+=degreeIncr;
    		((Vector2D)vertices[i]).E(x,y);
    	}
    	applyTranslationRotation();
    }
    
    private double edgeLength;
    private double n, p;//n = -size/2, p = +size/2
    private double radius;
    private double x,y;
    private int nsides;

}
