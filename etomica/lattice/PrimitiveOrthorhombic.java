package etomica.lattice;
import etomica.Space;

/**
 * Primitive group for an orthorhombic system.  All primitive
 * vectors orthogonal but not necessarily of equal length.
 */
public class PrimitiveOrthorhombic extends Primitive {
    
    public PrimitiveOrthorhombic(Space space) {
        super(space);
        //set up orthogonal vectors of unit size
        setSize(1.0);
    }
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.
     */
    public void setSize(double[] size) {
        if(size.length != D) throw new IllegalArgumentException("Error in PrimitiveOrthorhombic.setSize: Number of sizes given is inconsistent with number of primitive vectors");
        for(int i=0; i<D; i++) r[i].setComponent(i,size[i]);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        for(int i=0; i<D; i++) r[i].setComponent(i,size);
    }
    
    public Primitive reciprocal() {
        throw new RuntimeException("method PrimitiveOrthorhombic.reciprocal not yet implemented");
    }
    
    public AbstractCell wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public AbstractCell unitCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.unitCell not yet implemented");
    }
}//end of PrimitiveOrthorhombic
    
