/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.Space;
import etomica.space3d.Vector3D;

/**
 * A polyhedron composed of three pairs of rectangular faces placed opposite
 * each other and joined at right angles, also known as a
 * rectangular parallelepiped.
 * @author kofke
 * 
 */
public class Cuboid extends Hexahedron {

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
        return a * b * c;
    }

    public double getSurfaceArea() {
        return 2*(a*b + a*c + b*c);
    }
    
    public double getPerimeter() {
        return 4*(a + b + c);
    }

    public void updateVertices() {
        ((Vector3D)vertices[0]).E(na, nb, nc);
        ((Vector3D)vertices[1]).E(na, nb, pc);
        ((Vector3D)vertices[2]).E(na, pb, nc);
        ((Vector3D)vertices[3]).E(na, pb, pc);
        ((Vector3D)vertices[4]).E(pa, nb, nc);
        ((Vector3D)vertices[5]).E(pa, nb, pc);
        ((Vector3D)vertices[6]).E(pa, pb, nc);
        ((Vector3D)vertices[7]).E(pa, pb, pc);
        applyTranslationRotation();
    }

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the
     * surface of) this cell, <code>false</code> otherwise.
     */
    public boolean contains(etomica.space.Vector v) {
        double x = v.x(0)-position.x(0);
        double y = v.x(1)-position.x(1);
        double z = v.x(2)-position.x(2);
        return (x >= na) && (x <= pa) && (y >= nb) && (y <= pb) && (z >= nc)
                && (z <= pc);
    }
    
    /**
     * Sets the lengths of all edges of the cuboid.
     */
    public void setEdgeLengths(double a, double b, double c) {
        this.a = a;
        this.b = b;
        this.c = c;
        na = -0.5 * a;
        pa = +0.5 * a;
        nb = -0.5 * b;
        pb = +0.5 * b;
        nc = -0.5 * c;
        pc = +0.5 * c;
        updateVertices();
    }

    private double a, b, c;
    private double na, pa;//na = -a/2, p = +a/2
    private double nb, pb;//nb = -b/2, p = +a/2
    private double nc, pc;//nc = -c/2, p = +a/2

}
