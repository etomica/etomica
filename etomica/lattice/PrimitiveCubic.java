package etomica.lattice;
import etomica.Space;

public class PrimitiveCubic extends Primitive {
    
    public PrimitiveCubic(Space space) {
        super(space);
        //set up orthogonal vectors of unit size
        for(int i=0; i<D; i++) {
            r[i].setComponent(i,1.0);
        }
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
}
    
