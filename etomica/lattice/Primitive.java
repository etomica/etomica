package etomica.lattice;
import etomica.*;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
public abstract class Primitive {
    
    protected final Space.Vector[] latticeVectors;
    protected final Space.Vector[] latticeVectorsCopy;
    protected final int[] idx;//used to return coordinate index
    public final int D;
    protected Space space;
    protected Simulation simulation;
    protected BravaisLattice lattice;
    
    public Primitive(Simulation sim) {
        simulation = sim;
        space = sim.space;
        D = space.D();
        latticeVectors = new Space.Vector[D];
        latticeVectorsCopy = new Space.Vector[D];
        idx = new int[D];
        for(int i=0; i<D; i++) {
            latticeVectors[i] = space.makeVector();
            latticeVectorsCopy[i] = space.makeVector();
        }
    }
    
    /**
     * Lattice associated with this may be defined to enable
     * automatic updates of it with changes in the basis.
     * This method is called by the setPrimitive method of
     * BravaisLattice when the primitive is assigned to the lattice.
     */
    public void setLattice(BravaisLattice lattice) {
        this.lattice = lattice;
    }
    /**
     * Lattice associated with this primitive.
     */
    public BravaisLattice getLattice() {return lattice;}
    
    /**
     * Returns the primitive vectors.  Does not return the original
     * vectors used by the class, but instead returns copies.  Thus
     * changing the vectors returned by this method does not modify
     * the primitive vectors used by this class.  Subclasses should
     * provide mutator methods that permit changes to the vectors while
     * adhering to a particular structure (e.g., cubic, fcc, etc.).
     */
    public Space.Vector[] vectors() {
        return copyVectors();
    }
    
    //copies the interal set of vectors to the copy for outside use
    protected Space.Vector[] copyVectors() {
        for(int i=0; i<D; i++) latticeVectorsCopy[i].E(latticeVectors[i]);
        return latticeVectorsCopy;
    }
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public abstract Primitive copy();
    
    /**
     * Returns the index which would give the unit cell containing the given
     * point if the index were passed to a the site method of a sufficiently
     * large lattice that uses this primitive.
     */
    public abstract int[] latticeIndex(Space.Vector r);
    
    /**
     * Same as latticeIndex(Space.Vector), but gives index for periodic system
     * with number of unit cells in each direction as given by the dimensions array.
     * If lattice index corresponds to a cell outside the range of dimensions,
     * index of image in central cells is returned.
     */
    public abstract int[] latticeIndex(Space.Vector r, int[] dimensions);
    
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