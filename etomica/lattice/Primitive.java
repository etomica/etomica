package etomica.lattice;
import etomica.Space;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
public abstract class Primitive {
    
    protected final Space.Vector[] r;
    private final Space.Vector[] rCopy;
    public final int D;
    
    public Primitive(Space space) {
        D = space.D();
        r = new Space.Vector[D];
        rCopy = new Space.Vector[D];
        for(int i=0; i<D; i++) {
            r[i] = space.makeVector();
            rCopy[i] = space.makeVector();
        }
    }
    
    /**
     * Returns the primitive vectors.  Does not return the original
     * vectors used by the class, but instead returns copies.  Thus
     * changing the vectors returned by this method does not modify
     * the primitive vectors used by this class.  Subclasses should
     * provide mutator methods that permit changes to the vectors while
     * adhering to a particular structure (e.g., cubic, fcc, etc.).
     */
    public Space.Vector[] vectors() {
        return copy();
    }
    
    //copies the interal set of vectors to the copy for outside use
    private Space.Vector[] copy() {
        for(int i=0; i<D; i++) rCopy[i].E(r[i]);
        return rCopy;
    }
    
    /**
     * Returns the primitive for the reciprocal lattice vectors.
     */
    public abstract Primitive reciprocal();
    
    /**
     * Makes and returns a new Wigner-Seitz cell for the lattice 
     * specified by this primitive.
     */
    public abstract AbstractCell wignerSeitzCell();
    
    /**
     * Makes and returns a new unit cell for the lattice 
     * specified by this primitive.
     */
    public abstract AbstractCell unitCell();
    
}