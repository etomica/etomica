package etomica.lattice;

import etomica.math.geometry.Plane;
import etomica.*;

/**
 * Class describing a plane through a lattice.  Holds a Plane
 * object and a Primitive.  The orientation of the plane is 
 * specified using Miller indices defined according to the primitive.
 */
 
 /* History
  * 09/07/02 (DAK) new
  */
  
public class LatticePlane implements AtomFilter {
    
    private final Plane plane;
    private Primitive primitive;
    private Primitive reciprocal;
    private Space.Vector normal;
    private Space space;
    private int[] millerIndices;
    private int position;
    
    public LatticePlane(Primitive primitive, int[] h) {
        this.primitive = primitive;
        reciprocal = primitive.reciprocal();
        space = primitive.space;
        millerIndices = new int[space.D()];
        switch(space.D()) {
            case 3: plane = new Plane(); break;
            default: throw new IllegalArgumentException("LatticePlane not defined for other than 3D space");
        }
        normal = space.makeVector();
        setMillerIndices(h);
    }
    
    public boolean accept(Atom a) {
        return plane.isPositiveSide((Space3D.Vector)a.coord.position());
    }
    
    public void setPrimitive(Primitive primitive) {
        this.primitive = primitive;
        reciprocal = primitive.reciprocal();
        //update orientation of plane using current miller indices and the new primitive
        setMillerIndices(millerIndices);
    }
    
    /**
     * Sets the plane to the orientation specified by the given Miller indices
     * applied to the current primitive.
     */
    public void setMillerIndices(int[] h) {
        if(h.length != space.D()) throw new IllegalArgumentException("Error: number of miller indices passed to LatticePlane.setMillerIndices inconsistent with spatial dimension");
        normal.E(0.0);
        Space.Vector[] b = reciprocal.vectors();
        for(int i=0; i<h.length; i++) {
            normal.PEa1Tv1(h[i],b[i]);
            millerIndices[i] = h[i];
        }
        plane.setNormalVector(normal);
    }
    
    /**
     * Returns the current value of the i-th Miller index.
     */
    public int getMillerIndex(int i) {return millerIndices[i];}
    /**
     * Sets the i-th Miller index to h, keeping the others fixed.
     */
    public void setMillerIndex(int i, int h) {
        if(i < 0 || i >= space.D()) throw new IllegalArgumentException("Error: specified miller index passed to LatticePlane.setMillerIndex inconsistent with spatial dimension");
        millerIndices[i] = h;
        setMillerIndices(millerIndices);
    }
    
    /**
     * Sets to the position of the i-th plane from the one through the origin.
     */
    public void setPosition(int i) {
        position = i;
        //use normal as a work vector, giving a point contained by the plane
        //in its desired position
        Space.Vector[] b = reciprocal.vectors();
        for(int j=0; j<millerIndices.length; j++) {
            normal.PEa1Tv1(2.0*Math.PI*i*millerIndices[j],b[j]);
        }
        plane.moveTo((Space3D.Vector)normal);
    }
    
    /**
     * Returns the position i, such the the plane is located at the i-th from the origin.
     */
    public int getPosition() {return position;}

    /**
     * Returns true if the given point is on the side of the 
     * plane toward which the normal vector points.  The direction
     * of the normal vector can be inverted using the invert method.
     */
    public boolean isPositiveSide(Space3D.Vector p) {
        return plane.isPositiveSide(p);
    }
    
    /* implement this
    /**
     * Returns true if the given point is inside the plane (within some small tolerance).
     * /
    public boolean isInside(Space.Vector p) {
        return plane.isInside(p);
    }
    
    /**
     * Changes the direction of the normal vector so that it points
     * toward the other side of the plane from its present orientation.
     * Does not affect the location or absolute orientation of the plane.
     */
    public void invert() {
        plane.invert();
    }
}//end of LatticePlane