package simulate.lattice;
import simulate.Space;

/**
 * A Site with vertices and a volume.
 */
public abstract class AbstractCell extends Site {
    
    public AbstractCell(AbstractLattice parent, AbstractLattice.PositionCoordinate coord) {
        super(parent, new SiteIterator.Neighbor(), coord);
        //position null at this point
//        if(coord.position().D() != this.D()) { //define an exception for this (DimensionConflictException ?)
//            System.out.println("Dimension conflict in Cell constructor");
//            System.exit(1);
//        }
    }
    /**
     * Dimension of the space occupied by the cell
     */
     public abstract int D();
     
    /**
     * Returns the volume of the cell.
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
     * A vertex is a point where two or more edges meet.
     */
    public final int vertexCount() {return vertex().length;}
    
    /**
     * Returns squared distance between nearest vertices of this cell and the given cell.
     */
    public double r2NearestVertex(AbstractCell c, Space.Boundary boundary) {
        Space.Vector[] v1 = vertex();
        Space.Vector[] v2 = c.vertex();
        double r2Min = Double.MAX_VALUE;
        for(int k1=0; k1<vertexCount(); k1++) {
            for(int k2=0; k2<c.vertexCount(); k2++) {
                r2Min = Math.min(r2Min, Space.r2(v1[k1],v2[k2],boundary));
            }
        }
        return r2Min;
    }
}
