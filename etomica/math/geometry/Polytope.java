/*
 * History
 * Created on Nov 24, 2004 by kofke
 */
package etomica.math.geometry;

import etomica.Space;

/**
 * @author kofke
 *
 * Representation of a mathematical polytope, which is a finite region
 * of space enclosed by a finite number of hyperplanes.  Subclasses
 * include LineSegment (a 1-D polytope), Polygon (a 2-D polytope), and
 * Polyhedron (a 3-D polytope).
 */

//TODO add more general methods to specify polytope, and generate features such as vertices
public abstract class Polytope {

    public Polytope(Space space) {
        this.space = space;
    }
    
    /**
     * Specifies a linear size of the polytope.  For example, for 
     * a cube this is the length of each edge.  The definition of the
     * size will depend on the definition of the concrete subclass.
     * @param size
     */
    public abstract void setSize(double size);
    
    /**
     * Returns a specification of the linear size of the polytope.
     */
    public abstract double getSize();
     
    /**
     * Returns the (hyper)volume of the polytope.
     */
    public abstract double volume();
    
    /**
     * Returns the positions of the vertices relative to the cell position.
     * Absolute positions are obtained by adding the coordinate.position vector.
     * Note that vertices might be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values if using them often, 
     * but if doing so be careful to update them if any transformations are done to the lattice.
     */
    public abstract Space.Vector[] vertex();
    
    /**
     * Returns <code>true</code> if the given vector lies inside the cell, <code>false</code> otherwise.
     */
    public abstract boolean inCell(Space.Vector v);

    /**
     * Number of vertices bounding the cell.
     * A vertex is a point where D edges meet.
     */
    public final int vertexCount() {return vertex().length;}
    
    /**
     * Returns squared distance between nearest vertices of this cell
     * and an identical cell separated from it by the given vector.
     */
    public double r2NearestVertex(Space.Vector dr) {
        Space.Vector[] v = vertex();
        double r2Min = Double.MAX_VALUE;
        Space.Vector vdr = (Space.Vector)dr.clone();
        for(int k1=0; k1<v.length; k1++) {
            for(int k2=0; k2<v.length; k2++) {
                //compute [v[k1]-(v[k2]+dr)]^2
                vdr.Ev1Mv2(v[k1],v[k2]);               
                r2Min = Math.min(r2Min, vdr.Mv1Squared(dr));
                if(r2Min == 0.0) return r2Min;//minimum possible value -- no need to look at others
            }
        }
        return r2Min;
    }
    
    public Space space() {
        return space;
    }
    
    public final Space space;

}
