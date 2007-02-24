/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.exception.MethodNotImplementedException;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space2d.Vector2D;

/**
 * A quadrilateral with opposite sides parallel
 * .
 * @author kofke
 * 
 */
public class Parallelogram extends Polygon implements Parallelotope {

    /**
     * Default constructor makes a square of unit size
     */
    public Parallelogram(Space embeddedSpace) {
        this(embeddedSpace, new Vector2D(1,0), new Vector2D(0,1));
    }

    /**
     * Constructs a parallelogram with the given edge vectors.
     */
    public Parallelogram(Space embeddedSpace, Vector2D a, Vector2D b) {
        super(embeddedSpace, 4);
        this.a = (Vector2D)embeddedSpace.makeVector();
        this.b = (Vector2D)embeddedSpace.makeVector();
        setEdgeVectors(new IVector[] {a, b});
    }

    public double getArea() {
        return 4.0*Math.abs(a.x(0)*b.x(1) - a.x(1)*b.x(0));//4.0 is because a,b,c are half edge length
    }
    
    public double getPerimeter() {
        double aL = Math.sqrt(a.squared());
        double bL = Math.sqrt(b.squared());
        return 4 * (aL + bL);
    }

    /**
     * Calculated vertices based on current values of edge vectors.
     */
    public void updateVertices() {
        vertices[3].Ev1Mv2(a, b);//+a, -b
        vertices[2].Ev1Pv2(a, b);//+a, +b
        vertices[1].Ea1Tv1(-1, vertices[3]);//-a, +b
        vertices[0].Ea1Tv1(-1, vertices[2]);//-a, -b
        applyTranslationRotation();
    }

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the
     * surface of) this cell, <code>false</code> otherwise.
     */
    //TODO implement contains method in Parallelogram
    public boolean contains(IVector v) {
        throw new MethodNotImplementedException();
//        double x = v.x(0)-position.x(0);
//        double y = v.x(1)-position.x(1);
//        double z = v.x(2)-position.x(2);
//        return (x >= na) && (x <= pa) && (y >= nb) && (y <= pb) && (z >= nc)
//                && (z <= pc);
    }
    
    /**
     * Sets the lengths and directions of all edges of the parellelepiped.
     * Given instances are copied to an internal representation.
     */
    public void setEdgeVectors(IVector[] edgeVectors) {
        a.Ea1Tv1(0.5, edgeVectors[0]);
        b.Ea1Tv1(0.5, edgeVectors[1]);
        updateVertices();
    }

    private final Vector2D a, b;
    private static final long serialVersionUID = 1L;
}
