package etomica.lattice.crystal;
import etomica.math.geometry.Polytope;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Primitive group for a 2D hexagonal system.  Primitive-vector angle
 * is 120 degrees and two vectors are of equal length.
 */
public class PrimitiveHexagonal2D extends Primitive {
    
    private static final long serialVersionUID = 1L;
    protected double ab;
    protected static final double gamma = etomica.units.Degree.UNIT.toSim(120.);
    protected static final double cosGamma = Math.cos(gamma);
    protected static final double sinGamma = Math.sin(gamma);
    
    public PrimitiveHexagonal2D(Space space) {
        this(space, 1.0);
    }
    public PrimitiveHexagonal2D(Space space, double ab) {
        this(space, ab, true);
    }

    protected PrimitiveHexagonal2D(Space space, double ab, boolean makeReciprocal) {
        super(space, makeReciprocal);
        setSize(new double[]{ab, ab});
        setAngles(new double[]{gamma});
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveHexagonal2DReciprocal(space, 1);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        double[] newReciprocalSize = new double[2];
        newReciprocalSize[0] = 2.0 * Math.PI * size[0] * sinGamma;
        newReciprocalSize[1] = 2.0 * Math.PI * size[0] * sinGamma;
        reciprocal.setSize(newReciprocalSize);
    }
    
    public void setSizeAB(double newAB) {
        if (newAB == ab) {
            return;
        }
        setSize(new double[]{newAB, newAB});
        ab = newAB;
    }
    
    public double getSizeAB() {
        return ab;
    }
    
    protected void update() {
        super.update();
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(0,size[0]*cosGamma);
        latticeVectors[1].setX(1,size[0]*sinGamma);
    }
    
    /**
     * Returns a new PrimitiveHexagonal with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveHexagonal2D(space, ab);
    }
    
    public void scaleSize(double scale) {
        setSizeAB(scale*ab);
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
    
    public String toString() {return "Hexagonal";}
    

    protected static class PrimitiveHexagonal2DReciprocal extends PrimitiveHexagonal2D {
        public PrimitiveHexagonal2DReciprocal(Space space, double ab) {
            super(space, ab, false);
        }

        protected void update() {
            // this will screw up latticeVectors 0 and 1, but then we'll fix it
            super.update();
            latticeVectors[1].setX(1,ab);
            latticeVectors[0].setX(1,-ab*cosGamma);
            latticeVectors[0].setX(0,ab*sinGamma);
        }
        
        private static final long serialVersionUID = 1L;
    }
}
    
