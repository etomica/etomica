package etomica.lattice;
import etomica.Space;
import etomica.AtomFactory;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
public abstract class Primitive {
    
    protected final Space.Vector[] r;
    private final Space.Vector[] rCopy;
    protected final int[] idx;//used to return coordinate index
    public final int D;
    protected Space space;
    
    public Primitive(Space space) {
        this.space = space;
        D = space.D();
        r = new Space.Vector[D];
        rCopy = new Space.Vector[D];
        idx = new int[D];
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
     * Returns the index which would give the unit cell containing the given
     * point if the index were passed to a the site method of a sufficiently
     * large lattice that uses this primitive.
     */
    public abstract int[] latticeIndex(Space.Vector r);
    
    /**
     * Returns the primitive for the reciprocal lattice vectors.
     */
    public abstract Primitive reciprocal();
    
    /**
     * Returns a factory for the Wigner-Seitz cell specified by this primitive.
     * Each WS cell made by this factory remains tied to the primitive,
     * so its behavior will change to reflect any subsequent changes
     * in the primitive itself.
     */
    public abstract AtomFactory wignerSeitzCellFactory();
    
    /**
     * Returns a factory for the unit cell specified by this primitive.
     * Each unit cell made by this factory remains tied to the primitive,
     * so its behavior will change to reflect any subsequent changes
     * in the primitive itself.
     */
    public abstract AtomFactory unitCellFactory();
    
}