package etomica.lattice.crystal;
import etomica.lattice.Primitive;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Primitive group for a face-centered-cubic system.
 */
public class PrimitiveFcc extends Primitive implements Primitive3D {
    
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    //default size is 1.0
    private double size;
    private Vector[] unitVectors;
    private static final double FCC_ANGLE = Math.acos(0.5);
    
    public PrimitiveFcc(Space space) {
        this(space, 1.0);
    }
    public PrimitiveFcc(Space space, double size) {
        super(space); //also makes reciprocal
        //set up orthogonal vectors of unit size
        unitVectors = new Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
        setSize(size); //also sets reciprocal via update
    }
    /**
     * Constructor used by makeReciprocal method of PrimitiveBcc.
     */
    PrimitiveFcc(Space space, Primitive direct) {
        super(space, direct);
        unitVectors = new Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveBcc(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveBcc)reciprocal()).setSize(4.0*Math.PI/size);
    }
    
    public void setA(double a) {setSize(a);}
    public double getA() {return size;}
    
    public void setB(double b) {setSize(b);}
    public double getB() {return size;}
        
    public void setC(double c) {setSize(c);}
    public double getC() {return size;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return FCC_ANGLE;}
    
    public void setBeta(double t) {}
    public double getBeta() {return FCC_ANGLE;}
    
    public void setGamma(double t) {}
    public double getGamma() {return FCC_ANGLE;}
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return false;}
    public boolean isEditableC() {return false;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return false;}
    public boolean isEditableGamma() {return false;}

    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveFcc(space, size);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        this.size = size;
        for(int i=0; i<D; i++) latticeVectors[i].Ea1Tv1(size,unitVectors[i]);
        update();
    }
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getCubicSize() {return size;}
    
    /**
     * Multiplies the size of the current vectors by the given value.
     */
    public void scaleSize(double scale) {
        setSize(scale*size);
    }

    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
/*        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
*/    }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
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
    
    public Polytope unitCell() {
        throw new RuntimeException("method unitCellFactory not yet implemented");
    }
    
    public String toString() {return "Fcc";}

}
