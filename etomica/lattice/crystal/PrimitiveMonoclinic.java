package etomica.lattice.crystal;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Primitive group for a monoclinic system.  One primitive-vector
 * angle not normal, and vectors not necessarily of equal length.
 * a != b != c; alpha = gamma = 90deg, beta >= 90deg
 */
public class PrimitiveMonoclinic extends Primitive {
    
    private static final long serialVersionUID = 1L;

    public PrimitiveMonoclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle);
    }
    public PrimitiveMonoclinic(Space space, double a, double b, double c, double beta) {
        this(space, a, b, c, beta, true);
    }
    
    protected PrimitiveMonoclinic(Space space, double a, double b, double c, 
                               double beta, boolean makeReciprocal) {
        super(space, makeReciprocal);//also makes reciprocal
        setSize(new double[]{a, b, c});//also sets reciprocal via update
        setBeta(beta);
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveMonoclinic(space, 1, 1, 1, rightAngle, false);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveMonoclinic)reciprocal).setSize(new double[]{2.0*Math.PI/(size[0]*Math.sin(angle[1])),
                    2.0*Math.PI/size[1], 2.0*Math.PI/(size[2]*Math.sin(angle[1]))});
        ((PrimitiveMonoclinic)reciprocal).setBeta(angle[1]);
    }
    
    public void setA(double newA) {
        if (size[0] == newA) {
            return;
        }
        setSize(new double[]{newA, size[1], size[2]});
    }
    public double getA() {return size[0];}
    
    public void setB(double newB) {
        if (size[1] == newB) {
            return;
        }
        setSize(new double[]{size[0], newB, size[2]});
    }
    public double getB() {return size[1];}
        
    public void setC(double newC) {
        if (size[2] == newC) {
            return;
        }
        setSize(new double[]{size[0], size[1], newC});
    }
    public double getC() {return size[2];}

    //direct lattice (ix = 0, iz = 2)
    // v[0] = (1,0,0); v[1] = (0,1,0); v[2] = (c,0,s)  (times a, b, c)
     
    //reciprocal lattice (ix = 2, iz = 0)
    // v[0] = (s,0,-c); v[1] = (0,1,0); v[2] = (0,0,1);  (times a, b, c)
    protected void update() {
        super.update();
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(1,size[1]);
        latticeVectors[2].setX(0,size[2]*Math.cos(angle[1]));
        latticeVectors[2].setX(2,size[2]*Math.sin(angle[1]));
    }
    
    public void setBeta(double t) {
        if (t < rightAngle || t > Math.PI) {
            throw new IllegalArgumentException("Beta must be between PI/2 and PI");
        }
        setAngles(new double[]{rightAngle, t, rightAngle});
    }
    
    public double getBeta() {
        return angle[1];
    }

    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveMonoclinic(space, size[0], size[1], size[2], angle[1]);
    }
        
    public void scaleSize(double scale) {
        setSize(new double[]{size[0]*scale, size[1]*scale, size[2]*scale});
    }        
    
    public int[] latticeIndex(Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public Polytope unitCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public String toString() {return "Monoclinic";}

    protected static class PrimitiveMonoclinicReciprocal extends PrimitiveMonoclinic {
        private static final long serialVersionUID = 1L;

        public PrimitiveMonoclinicReciprocal(Space space, double a, double b, double c, double beta) {
            super(space, a, b, c, beta, false);
        }

        protected void updateLatticeVectors() {
            // this will screw up latticeVectors 0 and 1, but then we'll fix it
            super.update();
            latticeVectors[0].setX(0,size[0]*Math.sin(angle[1]));
            latticeVectors[0].setX(2,-size[0]*Math.cos(angle[1]));
            latticeVectors[1].setX(1,size[1]);
        }
    }
}
