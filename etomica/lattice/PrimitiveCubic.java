package etomica.lattice;
import etomica.Space;

/**
 * Primitive group for a cubic system.  All primitive
 * vectors orthogonal and of equal length.
 */
public class PrimitiveCubic extends Primitive {
    
    public PrimitiveCubic(Space space) {
        super(space);
        //set up orthogonal vectors of unit size
        setSize(1.0);
    }
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.
     */
/*    public void setSize(double[] size) {
        if(size.length != D) throw new IllegalArgumentException("Error in PrimitiveCubic.setSize: Number of sizes given is inconsistent with number of primitive vectors");
        for(int i=0; i<D; i++) r[i].setComponent(i,size[i]);
    }
*/    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        for(int i=0; i<D; i++) r[i].setComponent(i,size);
    }
    
    public Primitive reciprocal() {
        throw new RuntimeException("method PrimitiveCubic.reciprocal not yet implemented");
    }
    
    public AbstractCell wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveCubic.wignerSeitzCell not yet implemented");
    }
    
    public AbstractCell unitCell() {
        throw new RuntimeException("method PrimitiveCubic.unitCell not yet implemented");
    }
}//end of PrimitiveCubic
    
