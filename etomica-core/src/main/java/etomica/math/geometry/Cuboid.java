/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * A polyhedron composed of three pairs of rectangular faces placed opposite
 * each other and joined at right angles, also known as a
 * rectangular parallelepiped.
 * @author kofke
 * 
 */
public class Cuboid extends Hexahedron implements Rectangular {

    /**
     * Constructs a cuboid with equal faces of unit size (a cube).
     */
    public Cuboid(Space embeddedSpace) {
        this(embeddedSpace, 1.0, 1.0, 1.0);
    }

    /**
     * Constructs a cuboid with edges of lengths having the given values.
     */
    public Cuboid(Space embeddedSpace, double a, double b, double c) {
        super(embeddedSpace);
        setEdgeLengths(a, b, c);
    }

    public double getVolume() {
        return 8 * pa * pb * pc;
    }

    public double getSurfaceArea() {
        return 2*4*(pa*pb + pa*pc + pb*pc);
    }
    
    public double getPerimeter() {
        return 4*2*(pa + pb + pc);
    }

    public void updateVertices() {
        vertices[0].setX(0,-pa);
        vertices[0].setX(1,-pb);
        vertices[0].setX(2,-pc);
        vertices[1].E(vertices[0]);
        vertices[1].setX(2,pc);
        vertices[2].E(vertices[0]);
        vertices[2].setX(1,pb);
        vertices[3].E(vertices[2]);
        vertices[3].setX(2,pc);
        vertices[4].E(vertices[0]);
        vertices[4].setX(0,pa);
        vertices[5].E(vertices[4]);
        vertices[5].setX(2,pc);
        vertices[6].E(vertices[4]);
        vertices[6].setX(1,pb);
        vertices[7].E(vertices[6]);
        vertices[7].setX(2,pc);
        applyTranslationRotation();
    }

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the
     * surface of) this cell, <code>false</code> otherwise.
     */
    public boolean contains(Vector v) {
        double x = v.getX(0)-position.getX(0);
        double y = v.getX(1)-position.getX(1);
        double z = v.getX(2)-position.getX(2);
        return (x >= -pa) && (x <= pa) && (y >= -pb) && (y <= pb) && (z >= -pc)
                && (z <= pc);
    }
    
    /**
     * Sets the three edge lengths of the cuboid, taking
     * each element of the given vector for the length of
     * the corresponding cuboid edge.
     */
    public void setEdgeLengths(Vector e) {
        setEdgeLengths(e.getX(0), e.getX(1), e.getX(2));
    }
    
    /**
     * Returns a vector with elements equal to the
     * edge lengths of the cuboid.  The returned vector is not
     * used to represent the cuboid internally, so changing its
     * values will not affect the state of the cuboid.
     */
    public Vector getEdgeLengths() {
        return edgeLengths;
    }
    
    /**
     * Sets the lengths of all edges of the cuboid.
     */
    public void setEdgeLengths(double a, double b, double c) {
        edgeLengths.setX(0,a);
        edgeLengths.setX(1,b);
        edgeLengths.setX(2,c);
        pa = +0.5 * a;
        pb = +0.5 * b;
        pc = +0.5 * c;
        updateVertices();
    }

    private final Vector edgeLengths = new Vector3D();//used only to return edge lengths as a vector
    private double pa;//p = +a/2
    private double pb;//p = +a/2
    private double pc;//p = +a/2

}
