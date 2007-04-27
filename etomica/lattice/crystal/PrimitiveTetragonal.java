package etomica.lattice.crystal;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Polytope;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Primitive group for a tetragonal system.  All primitive
 * vectors orthogonal and two are of equal length.
 */
public class PrimitiveTetragonal extends Primitive {
    
    private static final long serialVersionUID = 1L;
    private double ab = 1.0;
    
    public PrimitiveTetragonal(Space space) {
        this(space, 1.0, 1.0);
    }
    public PrimitiveTetragonal(Space space, double ab, double c) {
        this(space, ab, c, true);
    }
    
    protected PrimitiveTetragonal(Space space, double ab, double c, boolean makeReciprocal) {
        super(space, makeReciprocal);//also makes reciprocal
        setSize(new double[]{ab, ab, c});
        setAngles(new double[]{rightAngle, rightAngle, rightAngle});
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveTetragonal(space, 1, 1, false);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        double[] newReciprocalSize = new double[3];
        newReciprocalSize[0] = 2.0 * Math.PI * size[0];
        newReciprocalSize[1] = 2.0 * Math.PI * size[0];
        newReciprocalSize[2] = 2.0 * Math.PI * size[2];
        reciprocal.setSize(newReciprocalSize);
    }
    
    public void setSizeAB(double newAB) {
        if (newAB == ab) {
            return;
        }
        setSize(new double[]{newAB, newAB, size[2]});
        ab = newAB;
    }

    public double getSizeAB() {
        return ab;
    }
    
    public void setSizeC(double newC) {
        if (newC == size[2]) {
            return;
        }
        setSize(new double[]{ab, ab, newC});
    }
    public double getSizeC() {return size[2];}

    public double getAngleAlpha() {
        return angle[0];
    }

    public double getAngleBeta() {
        return angle[1];
    }

    public double getAngleGamma() {
        return angle[2];
    }

    protected void update() {
        super.update();
        latticeVectors[0].setX(0, size[0]);
        latticeVectors[1].setX(1, size[1]);
        latticeVectors[2].setX(2, size[2]);
    }
    
    /**
     * Returns a new PrimitiveTetragonal with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveTetragonal(space, ab, size[2]);
    }
    
    public void scaleSize(double scale) {
        setSize(new double[]{ab*scale, ab*scale, size[2]*scale});
        ab = ab * scale;
    }

    public int[] latticeIndex(IVector q) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
   */ }
    
    public int[] latticeIndex(IVector q, int[] dimensions) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    */}
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method wignerSeitzCell not yet implemented");
    }
    
    /**
     * Returns a new Cuboid with edges given by the lengths of the
     * primitive vectors.
     */
    public Polytope unitCell() {
        return new Cuboid(space, ab, ab, size[2]);
    }
    
    public String toString() {return "Tetragonal";}
}
