/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

/**
 * A polyhedron composed of three pairs of rectangular faces placed opposite
 * each other and joined at right angles, also known as a
 * rectangular parallelepiped.
 * @author kofke
 * 
 */
public class Cuboid extends Polyhedron {

    /**
     * Constructs a cuboid with equal faces of unit size (a cube).
     */
    public Cuboid() {
        this(1.0, 1.0, 1.0);
    }

    /**
     * Constructs a cube with edges of lengths having the given values.
     */
    public Cuboid(double a, double b, double c) {
        super();
        vertices = new etomica.space3d.Vector[8];
        for (int i = 0; i < vertices.length; i++)
            vertices[i] = new etomica.space3d.Vector();
        setSize(a, b, c);
    }

    /**
     * Returns size^3.
     */
    public double volume() {
        return a * b * c;
    }

    /**
     * Returns the absolute positions of the vertices . Note that vertices might
     * be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values
     * if using them often, but if doing so be careful to update them if any
     * transformations are done to the lattice.
     */
    public etomica.space.Vector[] vertex() {
        vertices[0].E(na, nb, nc);
        vertices[1].E(na, nb, pc);
        vertices[2].E(na, pb, nc);
        vertices[3].E(na, pb, pc);
        vertices[4].E(pa, nb, nc);
        vertices[5].E(pa, nb, pc);
        vertices[6].E(pa, pb, nc);
        vertices[7].E(pa, pb, pc);
        return vertices;
    }//end of vertex

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the
     * surface of) this cell, <code>false</code> otherwise.
     */
    public boolean inCell(etomica.space.Vector v) {
        double x = v.x(0);
        double y = v.x(1);
        double z = v.x(2);
        return (x >= na) && (x <= pa) && (y >= nb) && (y <= pb) && (z >= nc)
                && (z <= pc);
    }

    /**
     * @return Returns the size, which is the length the first edge of the cuboid.
     */
    public double getSize() {
        return a;
    }

    /**
     * Sets the value of the first length, and scales the others
     * to maintain the shape of the cuboid.
     */
    public void setSize(double size) {
        double scale = size/a;
        setSize(a, scale*b, scale*c);
    }
    
    /**
     * Sets the lengths of all edges of the cuboid.
     */
    public void setSize(double a, double b, double c) {
        this.a = a;
        this.b = b;
        this.c = c;
        na = -0.5 * a;
        pa = +0.5 * a;
        nb = -0.5 * b;
        pb = +0.5 * b;
        nc = -0.5 * c;
        pc = +0.5 * c;
    }

    private double a, b, c;
    private double na, pa;//na = -a/2, p = +a/2
    private double nb, pb;//nb = -b/2, p = +a/2
    private double nc, pc;//nc = -c/2, p = +a/2

    private final etomica.space3d.Vector[] vertices;

}
