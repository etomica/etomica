/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.space.Vector;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;

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
        this(embeddedSpace, Vector.of(new double[]{1, 0}), Vector.of(new double[]{0, 1}));
    }

    /**
     * Constructs a parallelogram with the given edge vectors.
     */
    public Parallelogram(Space embeddedSpace, Vector a, Vector b) {
        super(embeddedSpace, 4);
        this.a = embeddedSpace.makeVector();
        this.b = embeddedSpace.makeVector();
        setEdgeVectors(new Vector[] {a, b});
    }

    public double getArea() {
        return 4.0*Math.abs(a.getX(0)*b.getX(1) - a.getX(1)*b.getX(0));//4.0 is because a,b,c are half edge length
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
    public boolean contains(Vector v) {
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
    public void setEdgeVectors(Vector[] edgeVectors) {
        a.Ea1Tv1(0.5, edgeVectors[0]);
        b.Ea1Tv1(0.5, edgeVectors[1]);
        updateVertices();
    }

    private final Vector a, b;
    private static final long serialVersionUID = 1L;
}
