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
import etomica.space3d.Vector3D;

/**
 * A polyhedron composed of three pairs of rectangular faces placed opposite
 * each other and joined at arbitrary angles.
 * @author kofke
 * 
 */
public class Parallelepiped extends Hexahedron implements Parallelotope {

    /**
     * Default constructor makes a cube of unit size
     */
    public Parallelepiped(Space embeddedSpace) {
        this(embeddedSpace, new Vector3D(1,0,0), new Vector3D(0,1,0), new Vector3D(0,0,1));
    }

    /**
     * Constructs a parallelepiped with the given edge vectors.
     */
    public Parallelepiped(Space embeddedSpace, Vector a, Vector b, Vector c) {
        super(embeddedSpace);
        this.a = embeddedSpace.makeVector();
        this.b = embeddedSpace.makeVector();
        this.c = embeddedSpace.makeVector();
        work = embeddedSpace.makeVector();
        setEdgeVectors(new Vector[] {a, b, c});
    }

    /**
     * Calculated parallelepiped volume using (a X b) . c
     */
    public double getVolume() {
        work.E(a);
        work.XE(b);
        return 8.0*Math.abs(work.dot(c));//8.0 is because a,b,c are half edge length
    }

    public double getSurfaceArea() {
        double a2 = a.squared();
        double b2 = b.squared();
        double c2 = c.squared();
        double ab = a.dot(b);
        double ac = a.dot(c);
        double bc = b.dot(c);
        return 2*4*( Math.sqrt(a2*b2-ab*ab) +  //4.0 is becase a,b,c are half edge lengths
                   Math.sqrt(a2*c2-ac*ac) + 
                   Math.sqrt(b2*c2-bc*bc));
    }
    
    public double getPerimeter() {
        double aL = Math.sqrt(a.squared());
        double bL = Math.sqrt(b.squared());
        double cL = Math.sqrt(c.squared());
        return 4 * (aL + bL + cL);
    }

    /**
     * Calculated vertices based on current values of edge vectors.
     */
    public void updateVertices() {
        vertices[7].Ev1Pv2(a, b);
        vertices[6].E(vertices[7]);
        vertices[5].Ev1Mv2(a, b);
        vertices[4].E(vertices[5]);
        vertices[7].PE(c);// +a, +b, +c
        vertices[6].ME(c);// +a, +b, -c
        vertices[5].PE(c);// +a, -b, +c
        vertices[4].ME(c);// +a, -b, -c
        vertices[3].Ea1Tv1(-1,vertices[4]);//-a, +b, +c
        vertices[2].Ea1Tv1(-1,vertices[5]);//-a, +b, -c
        vertices[1].Ea1Tv1(-1,vertices[6]);//-a, -b, +c
        vertices[0].Ea1Tv1(-1,vertices[7]);//-a, -b, -c
        applyTranslationRotation();
    }

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the
     * surface of) this cell, <code>false</code> otherwise.
     */
    //TODO implement contains method in Parallelepiped
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
    public void setEdgeVectors(Vector[] vectors) {
        a.Ea1Tv1(0.5, vectors[0]);
        b.Ea1Tv1(0.5, vectors[1]);
        c.Ea1Tv1(0.5, vectors[2]);
        updateVertices();
    }

    private static final long serialVersionUID = 1L;
    private final Vector a, b, c;
    private final Vector work;

}
