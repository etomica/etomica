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
    protected final double[] size;
    protected final double[] angle;
    private final double[] sizeCopy;
    protected Space space;
    protected Simulation simulation;
    protected BravaisLattice lattice;
    protected static final double rightAngle = 0.5*Math.PI;
    
    public Primitive(Simulation sim) {
        simulation = sim;
        space = sim.space;
        D = space.D();
        if(!( (this instanceof Primitive2D && D==2) || (this instanceof Primitive3D && D==3))) throw new RuntimeException("Error: inconsistency between spatial dimension and interface of Primitive");
        latticeVectors = new Space.Vector[D];
        latticeVectorsCopy = new Space.Vector[D];
        idx = new int[D];
        size = new double[D];
        sizeCopy = new double[D];
        angle = new double[D];
        for(int i=0; i<D; i++) {
            latticeVectors[i] = space.makeVector();
            latticeVectorsCopy[i] = space.makeVector();
            angle[i] = rightAngle;
        }
    }

    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.  Calls set[ABC] methods (defined in subclass)
     * for any lengths that are not equal to current values.
     */
    public void setSize(double[] newSize) {
        double size0, size1, size2;
        switch(D) {
            case 2:
                Primitive2D p2 = (Primitive2D)this;
                size0 = size[0];//save because might change with any of the set calls
                size1 = size[1];
                if(size0 != newSize[0]) p2.setA(newSize[0]);
                if(size1 != newSize[1]) p2.setB(newSize[1]);
                break;
            case 3:
                Primitive3D p3 = (Primitive3D)this;
                size0 = size[0];//save because might change with any of the set calls
                size1 = size[1];
                size2 = size[2];
                if(size0 != newSize[0]) p3.setA(newSize[0]);
                if(size1 != newSize[1]) p3.setB(newSize[1]);
                if(size2 != newSize[2]) p3.setC(newSize[2]);
                break;
            default:
                throw new RuntimeException("Didn't expect to get here in Primitive.setSize");
        }
        if(lattice != null) lattice.update();
    }
    
    /**
     * Sets the angles between the primitive vector to the corresponding
     * values in the given array.  Calls set[alpha/beta/gamma] methods (defined in subclass)
     * for any angles that are not equal to current values.
     */
    public void setAngles(double[] newAngle) {
        double t0, t1, t2;
        switch(D) {
            case 2:
                Primitive2D p2 = (Primitive2D)this;
                if(angle[0] != newAngle[0]) p2.setAlpha(newAngle[0]);
                break;
            case 3:
                Primitive3D p3 = (Primitive3D)this;
                t0 = angle[0];//save because might change with any of the set calls
                t1 = angle[1];
                t2 = angle[2];
                if(t0 != newAngle[0]) p3.setAlpha(newAngle[0]);
                if(t1 != newAngle[1]) p3.setBeta(newAngle[1]);
                if(t2 != newAngle[2]) p3.setGamma(newAngle[2]);
                break;
            default:
                throw new RuntimeException("Didn't expect to get here in Primitive.setSize");
        }
        if(lattice != null) lattice.update();
    }

    /**
     * Returns a copy of the array of primitive-vector sizes.
     */
/*    public double[] getSize() {
        for(int i=0; i<D; i++) sizeCopy[i] = size[i];
        return sizeCopy;
    }
*/    
    
    /**
     * Scales (multiplies) the size of each primitive vector by the given value.
     */
    public abstract void scaleSize(double scale);
    
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