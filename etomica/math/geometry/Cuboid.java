/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.api.IVector;
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
    public boolean contains(IVector v) {
        double x = v.x(0)-position.x(0);
        double y = v.x(1)-position.x(1);
        double z = v.x(2)-position.x(2);
        return (x >= na) && (x <= pa) && (y >= nb) && (y <= pb) && (z >= nc)
                && (z <= pc);
    }
    
    /**
     * Sets the three edge lengths of the cuboid, taking
     * each element of the given vector for the length of
     * the corresponding cuboid edge.
     */
    public void setEdgeLengths(IVector e) {
        setEdgeLengths(e.x(0), e.x(1), e.x(2));
    }
    
    /**
     * Returns a vector with elements equal to the
     * edge lengths of the cuboid.  The returned vector is not
     * used to represent the cuboid internally, so changing its
     * values will not affect the state of the cuboid.
     */
    public IVector getEdgeLengths() {
        return edgeLengths;
    }
    
    /**
     * Sets the lengths of all edges of the cuboid.
     */
    public void setEdgeLengths(double a, double b, double c) {
        edgeLengths.E(a, b, c);
        na = -0.5 * a;
        pa = +0.5 * a;
        nb = -0.5 * b;
        pb = +0.5 * b;
        nc = -0.5 * c;
        pc = +0.5 * c;
        updateVertices();
    }

    private final Vector3D edgeLengths = new Vector3D();//used only to return edge lengths as a vector
    private double na, pa;//na = -a/2, p = +a/2
    private double nb, pb;//nb = -b/2, p = +a/2
    private double nc, pc;//nc = -c/2, p = +a/2

}
