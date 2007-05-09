package etomica.lattice.crystal;
import etomica.math.geometry.Polytope;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Primitive group for a face-centered-cubic system.
 */
public class PrimitiveFcc extends Primitive {
    
    private static final long serialVersionUID = 1L;
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    //default size is 1.0
    private double cubicSize;
    private IVector[] unitVectors;
    private static final double FCC_ANGLE = Math.acos(0.5);
    
    public PrimitiveFcc(Space space) {
        this(space, 1.0);
    }
    public PrimitiveFcc(Space space, double size) {
        this(space, size, true);
    }
    
    protected PrimitiveFcc(Space space, double size, boolean makeReciprocal) {
        super(space, makeReciprocal); //also makes reciprocal
        //set up orthogonal vectors of unit size
        unitVectors = new IVector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
        setCubicSize(size); //also sets reciprocal via update
        double[] newAngles = new double[D];
        for (int i=0; i<D; i++) {
            newAngles[i] = FCC_ANGLE;
        }
        setAngles(newAngles);
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveBcc(space, 1, false);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveBcc)reciprocal).setCubicSize(Math.sqrt(6)*Math.PI/size[0]);
    }
    
    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveFcc(space, cubicSize);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setCubicSize(double newCubicSize) {
        if (newCubicSize == cubicSize) {
            // no change
            return;
        }
        double[] sizeArray = new double[D];
        for(int i=0; i<D; i++) {
            sizeArray[i] = newCubicSize;
        }
        setSize(sizeArray);
        cubicSize = newCubicSize;
    }

    protected void update() {
        super.update();
        for(int i=0; i<D; i++) latticeVectors[i].Ea1Tv1(size[0],unitVectors[i]);
    }
    
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getCubicSize() {return cubicSize;}
    
    /**
     * Multiplies the size of the current vectors by the given value.
     */
    public void scaleSize(double scale) {
        setCubicSize(scale*cubicSize);
    }

    public int[] latticeIndex(IVector q) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
/*        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
*/    }
    
    public int[] latticeIndex(IVector q, int[] dimensions) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
 /*       for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
 */   }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveFcc.wignerSeitzCell not yet implemented");
    }
    
    public String toString() {return "Fcc";}

}
