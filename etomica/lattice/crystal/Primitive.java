package etomica.lattice.crystal;

import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Parallelogram;
import etomica.math.geometry.Polytope;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.IVector3D;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
public abstract class Primitive implements java.io.Serializable {
    
    protected final IVector[] latticeVectors;
    protected final IVector[] latticeVectorsCopy;
    protected final int[] idx;//used to return coordinate index
    protected final int D;
    protected final double[] size;
    protected final double[] angle;
    protected final Space space;
    protected static final double rightAngle = 0.5*Math.PI;
    protected final Primitive reciprocal;
    
    /**
     * This constructor is called directly when a Primitive is constructing
     * its reciprocal primitive.  For construction of the direct-lattice
     * primitive, this constructor is called via the Primitive(Simulation) constructor.
     */
    protected Primitive(Space space, boolean makeReciprocal) {
        this.space = space;
        D = space.D();
        latticeVectors = new IVector[D];
        latticeVectorsCopy = new IVector[D];
        idx = new int[D];
        size = new double[D];
//        sizeCopy = new double[D];
        angle = new double[D];
        for(int i=0; i<D; i++) {
            latticeVectors[i] = space.makeVector();
            latticeVectorsCopy[i] = space.makeVector();
            angle[i] = rightAngle;
        }
        //if reciprocal is not null, this is a direct primitive; if it is null,
        //this is a reciprocal primitive of another primitive that is in the
        //process of being constructed.
        this.reciprocal = (makeReciprocal  ? makeReciprocal() : null);
    }
    
    /**
     * Method defining and constructing reciprocal primitive; called by
     * constructor of Primitive.  
     */
     //definition of this method should take care not to lead to calling of update 
     //method of the reciprocal primitive.
    protected abstract Primitive makeReciprocal();
    
    /**
     * Updates reciprocal primitive so that it is consistent with the
     * current parameters of this primitive.  Called by update method.
     */
    protected abstract void updateReciprocal();

    /**
     * @return the space
     */
    public final Space getSpace() {
        return space;
    }

    /**
     * Returns a new array with values equal to the lengths of the primitive vectors.
     */
    public double[] getSize() {
        return (double[])size.clone();
    }
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.
     */
    protected final void setSize(double[] newSize) {
        for (int i=0; i<newSize.length; i++) {
            if (newSize[i] <= 0.0) { 
                throw new IllegalArgumentException("sizes must be positive");
            }
        }
        for (int i=0; i<newSize.length; i++) {
            size[i] = newSize[i];
        }
        update();
    }
    
    /**
     * Sets the angles between the primitive vector to the corresponding
     * values in the given array.
     */
    protected final void setAngles(double[] newAngle) {
        for (int i=0; i<newAngle.length; i++) {
            if (newAngle[i] < 0 || newAngle[i] > Math.PI) {
                throw new IllegalArgumentException("Angles must be between 0 and pi");
            }
        }
        for (int i=0; i<newAngle.length; i++) {
            angle[i] = newAngle[i];
        }
        update();
    }

    protected void update() {
        for (int i=0; i<D; i++) {
            if (size[i] == 0) {
                // we haven't been fully set up yet
                return;
            }
        }
        for (int i=0; i<angle.length; i++) {
            if (angle[i] == 0) {
                // we haven't been fully set up yet
                return;
            }
        }
        if (reciprocal != null) {
            updateReciprocal();
        }
    }
        
    /**
     * Scales (multiplies) the size of each primitive vector by the given value.
     */
    public abstract void scaleSize(double scale);
        
    /**
     * Returns the primitive vectors.  Does not return the original
     * vectors used by the class, but instead returns copies.  Thus
     * changing the vectors returned by this method does not modify
     * the primitive vectors used by this class.  Subclasses should
     * provide mutator methods that permit changes to the vectors while
     * adhering to a particular structure (e.g., cubic, fcc, etc.).
     */
    public IVector[] vectors() {
        return copyVectors();
    }
    
    //copies the interal set of vectors to the copy for outside use
    protected IVector[] copyVectors() {
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
    public abstract int[] latticeIndex(IVector r);
    
    /**
     * Same as latticeIndex(Space.Vector), but gives index for periodic system
     * with number of unit cells in each direction as given by the dimensions array.
     * If lattice index corresponds to a cell outside the range of dimensions,
     * index of image in central cells is returned.
     */
    public abstract int[] latticeIndex(IVector r, int[] dimensions);
    
    /**
     * Returns the primitive for the reciprocal lattice vectors.
     */
    public Primitive reciprocal() {
        return reciprocal;
     }
        
    /**
     * Returns the Wigner-Seitz cell specified by this primitive.
     * The returned cell does not remain tied to the primitive, and
     * will not be updated with changes to the primitive.
     */
    public abstract Polytope wignerSeitzCell();
    
    /**
     * Returns a the unit cell specified by this primitive.
     * The returned cell does not remain tied to the primitive, and
     * will not be updated with changes to the primitive.
     */
    public Polytope unitCell() {
        if (space.D() == 1) {
            LineSegment line = new LineSegment(space);
            line.setLength(latticeVectors[0].x(0));
            return line;
        }
        if (space.D() == 2) {
            return new Parallelogram(space, latticeVectors[0], latticeVectors[1]);
        }
        if (space.D() == 3) {
            return new Parallelepiped(space, (IVector3D)latticeVectors[0], (IVector3D)latticeVectors[1], (IVector3D)latticeVectors[2]);
        }
        throw new RuntimeException("I'm impressed by your ability to make a "+D+"-D space, but I really don't know how to make an appropriate unit cell");
    }
    
    
}
