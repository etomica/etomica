/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

/**
 * A geometric cube.
 * 
 * @author kofke
 *  
 */
public class Cube extends Polyhedron {

    /**
     * Constructs a cube of unit size.
     */
    public Cube() {
        this(1.0);
    }
    
    /**
     * Constructs a cube with edge length having the given value.
     * @param size edge length of the cube
     */
    public Cube(double size) {
        super();
        vertices = new etomica.space3d.Vector[8];
        for(int i=0; i<vertices.length; i++) vertices[i] = new etomica.space3d.Vector();
        setSize(size);
    }

    /**
     * Returns size^3.
     */
    public double volume() {
        return size*size*size;
    }

    /**
     * Returns the absolute positions of the vertices .
     * Note that vertices might be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values if using them often, 
     * but if doing so be careful to update them if any transformations are done to the lattice.
     */
    public etomica.space.Vector[] vertex() {
        vertices[0].E(n,n,n);
        vertices[1].E(n,n,p);
        vertices[2].E(n,p,n);
        vertices[3].E(n,p,p);
        vertices[4].E(p,n,n);
        vertices[5].E(p,n,p);
        vertices[6].E(p,p,n);
        vertices[7].E(p,p,p);
        return vertices;
    }//end of vertex

    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean inCell(etomica.space.Vector v) {
        double x = v.x(0);
        double y = v.x(1);
        double z = v.x(2);
        return (x>=n) && (x<=p) && (y>=n) && (y<=p) && (z>=n) && (z<=p);
    }
    
    /**
     * @return Returns the size, which is the length the edge of the cube.
     */
    public double getSize() {
        return size;
    }
    /**
     * @param size The size to set.
     */
    public void setSize(double size) {
        this.size = size;
        n = -0.5*size;
        p = +0.5*size;
    }
    private double size;
    private double n, p;//n = -size/2, p = +size/2
    private final etomica.space3d.Vector[] vertices;

}
