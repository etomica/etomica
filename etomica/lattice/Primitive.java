package etomica.lattice;
import etomica.Space;

public abstract class Primitive {
    
    public final Space.Vector[] r;
    public final int D;
    
    public Primitive(Space space) {
        D = space.D();
        r = new Space.Vector[D];
        for(int i=0; i<D; i++) {r[i] = space.makeVector();}
    }
    
    public abstract Primitive reciprocal();
    
    public abstract AbstractCell wignerSeitzCell();
    
    public abstract AbstractCell unitCell();
    
}