package etomica.lattice;

import etomica.Space;
import etomica.math.geometry.Polytope;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
//TODO add propertyChangeListeners to be notified if primitive changes
public abstract class Primitive {
    
    protected final Space.Vector[] latticeVectors;
    protected final Space.Vector[] latticeVectorsCopy;
    protected final int[] idx;//used to return coordinate index
    public final int D;
    protected final double[] size;
    protected final double[] angle;
//    private final double[] sizeCopy;
    public final Space space;
    protected boolean immutable = false;//flag used when sync-ing with the reciprocal
    protected static final double rightAngle = 0.5*Math.PI;
    private final Primitive reciprocal;
    
    /**
     * This constructor is used by when making the direct-lattice primitive.
     */
    protected Primitive(Space space) {
        this(space, null);
    }
    /**
     * This constructor is called directly when a Primitive is constructing
     * its reciprocal primitive.  For construction of the direct-lattice
     * primitive, this constructor is called via the Primitive(Simulation) constructor.
     */
    protected Primitive(Space space, Primitive reciprocal) {
        this.space = space;
        D = space.D();
        if(!( (this instanceof Primitive2D && D==2) || (this instanceof Primitive3D && D==3))) throw new RuntimeException("Error: inconsistency between spatial dimension and interface of Primitive");
        latticeVectors = new Space.Vector[D];
        latticeVectorsCopy = new Space.Vector[D];
        idx = new int[D];
        size = new double[D];
//        sizeCopy = new double[D];
        angle = new double[D];
        for(int i=0; i<D; i++) {
            latticeVectors[i] = space.makeVector();
            latticeVectorsCopy[i] = space.makeVector();
            angle[i] = rightAngle;
        }
        //if reciprocal is null, this is a direct primitive; if it is not null,
        //this is a reciprocal primitive of another primitive that is in the
        //process of being constructed.
        this.reciprocal = (reciprocal == null) ? makeReciprocal() : reciprocal;
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
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.  Calls set[ABC] methods (defined in subclass)
     * for any lengths that are not equal to current values.
     */
    public void setSize(double[] newSize) {
        if(immutable) return;
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
        update();
    }
    
    /**
     * Sets the angles between the primitive vector to the corresponding
     * values in the given array.  Calls set[alpha/beta/gamma] methods (defined in subclass)
     * for any angles that are not equal to current values.
     */
    public void setAngles(double[] newAngle) {
        if(immutable) return;
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
        update();
    }

    /**
     * Returns a copy of the array of primitive-vector sizes.
     */
    //deleted because getSize is defined to return double in PrimitiveBcc, etc.
//    public double[] getSize() {
//        for(int i=0; i<D; i++) sizeCopy[i] = size[i];
//        return sizeCopy;
//    }
    

    protected void update() {
        if(immutable) return;
//        if(lattice != null) lattice.update();
        immutable = true; 
        updateReciprocal();
        immutable = false;
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
    public Primitive reciprocal() {return reciprocal;}
        
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
    public abstract Polytope unitCell();
    
}