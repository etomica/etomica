package etomica.lattice;
import etomica.*;

/**
 * Primitive group for an orthorhombic system.  All primitive
 * vectors orthogonal but not necessarily of equal length.
 */
public class PrimitiveOrthorhombic extends Primitive {
    
    public PrimitiveOrthorhombic(Simulation sim) {
        super(sim);
        size = new double[sim.space.D()];
        //set up orthogonal vectors of unit size
        setSize(1.0);
    }
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        PrimitiveOrthorhombic copy = new PrimitiveOrthorhombic(simulation);
        copy.setSize(size);
        return copy;
    }
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.
     */
    public void setSize(double[] size) {
        if(size.length != D) throw new IllegalArgumentException("Error in PrimitiveOrthorhombic.setSize: Number of sizes given is inconsistent with number of primitive vectors");
        for(int i=0; i<D; i++) {
            r[i].setComponent(i,size[i]);
            this.size[i] = size[i];
        }
        if(lattice != null) lattice.update();
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        for(int i=0; i<D; i++) {
            r[i].setComponent(i,size);
            this.size[i] = size;
        }
        if(lattice != null) lattice.update();
    }
    
    public int[] latticeIndex(Space.Vector q) {
        for(int i=0; i<D; i++) idx[i] = (int)(q.component(i)/size[i]);
        return idx;
    }

    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            idx[i] = (int)(q.component(i)/size[i]);
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Primitive reciprocal() {
        throw new RuntimeException("method PrimitiveOrthorhombic.reciprocal not yet implemented");
    }
    
    public AtomFactory wignerSeitzCellFactory() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public AtomFactory unitCellFactory() {
        throw new RuntimeException("method PrimitiveOrthorhombic.unitCell not yet implemented");
    }
    
    private double[] size;

///////////////////////////////////////////////////////////////////////////////////////////

public class UnitCellFactory extends AtomFactory {

    AtomType atomType;
    
    public UnitCellFactory(Simulation sim) {
        super(sim);
        setType(new AtomType(this));//default
    }
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}

    /**
     * Builds a single unit cell.
     */
    protected Atom build(AtomTreeNodeGroup parent) {
        return new UnitCell(space, atomType, parent);
    }

    public Atom build(Atom atom) {
        if(!(atom instanceof UnitCell)) throw new IllegalArgumentException("PrimitiveCubic.UnitCellFactory.build(Atom) attempted using an atom that is not an instance of UnitCell");
        return atom;
    }
    
}//end of UnitCellFactory

///////////////////////////////////////////////////////////////////////////////////////////

/**
 * A cubic unit cell.  Position of the cell is given by the vertex
 * in which each coordinate is minimized.
 */
public class UnitCell extends AbstractCell {
    
    public UnitCell(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, parent);
    }
    /**
     * Dimension of the space occupied by the cell
     */
     public int D() {return space.D();}
     
    /**
     * Returns the volume of the cubic cell.
     */
    public double volume() {
 /*       double sizeN = size;
        for(int i=D()-1; i>0; i--) sizeN *= size;
        return sizeN;*/
        return 0.0;
    }
    /**
     * Returns the positions of the vertices relative to the cell position.
     * Absolute positions are obtained by adding the coordinate.position vector.
     * Note that vertices might be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values if using them often, 
     * but if doing so be careful to update them if any transformations are done to the lattice.
     */
    public Space.Vector[] vertex() {
        return null;
    }
    
    /**
     * Returns <code>true</code> if the given vector lies inside the cell, <code>false</code> otherwise.
     */
    public boolean inCell(Space.Vector v) {
  /*      double x = size;
        switch(D()) {
            case 3: x = v.component(2);
                    if(x < 0.0 || x > r[2].component(2)) return false;
            case 2: x = v.component(1);
                    if(x < 0.0 || x > r[1].component(1)) return false;
            case 1: x = v.component(0);
                    if(x < 0.0 || x > r[0].component(0)) return false;
                    break;
            default: throw new RuntimeException("PrimitiveCubic.UnitCell.inCell not implemented for given dimension");
        }*/
        return true;
    }
}//end of UnitCell

}//end of PrimitiveOrthorhombic
    
