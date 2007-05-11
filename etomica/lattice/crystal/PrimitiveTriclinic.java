package etomica.lattice.crystal;
import etomica.math.geometry.Polytope;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Primitive group for a triclinic system.  No restrictions on
 * primitive-vector angles or lengths.
 * a != b != c; alpha != gamma != beta
 */
public class PrimitiveTriclinic extends Primitive {
    
    private static final long serialVersionUID = 1L;

    public PrimitiveTriclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle, rightAngle, rightAngle);
    }
    public PrimitiveTriclinic(Space space, double a, double b, double c, 
                                              double alpha, double beta, double gamma) {
        super(space);
        setSize(new double[]{a, b, c});
        setAngles(new double[]{alpha, beta, gamma});
    }

    //called by superclass constructor
    public Primitive makeReciprocal() {
        //XXX this does not update the reciprocal's size
//        PrimitiveTriclinic recip = (PrimitiveTriclinic)reciprocal;
//        Vector3D aStar = (Vector3D)recip.latticeVectors[0];
//        Vector3D bStar = (Vector3D)recip.latticeVectors[1];
//        Vector3D cStar = (Vector3D)recip.latticeVectors[2];
//        Vector3D aVec = (Vector3D)latticeVectors[0];
//        Vector3D bVec = (Vector3D)latticeVectors[1];
//        Vector3D cVec = (Vector3D)latticeVectors[2];
//        aStar.E(bVec);
//        aStar.XE(cVec);
//        double factor = 2.0*Math.PI/aVec.dot(aStar); // a . (b X c)
//        aStar.TE(factor);
//        bStar.E(cVec);
//        bStar.XE(aVec);
//        bStar.TE(factor);
//        cStar.E(aVec);
//        cStar.XE(bVec);
//        cStar.TE(factor);
        throw new RuntimeException("I don't know how to actually make the reciprocal for you.");
    }
    
    public void setSizeA(double newA) {
        if (size[0] == newA) {
            return;
        }
        setSize(new double[]{newA, size[1], size[2]});
    }
    public double getSizeA() {return size[0];}
    
    public void setSizeB(double newB) {
        if (size[1] == newB) {
            return;
        }
        setSize(new double[]{size[0], newB, size[2]});
    }
    public double getSizeB() {return size[1];}
        
    public void setSizeC(double newC) {
        if (size[2] == newC) {
            return;
        }
        setSize(new double[]{size[0], size[1], newC});
    }
    public double getSizeC() {return size[2];}

    protected void update() {
        super.update();
        double cosAlpha = Math.cos(angle[0]);
        double cosBeta = Math.cos(angle[1]);
        double cosGamma = Math.cos(angle[2]);
        double sinGamma = Math.sin(angle[2]);
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(0,size[1]*cosGamma);
        latticeVectors[1].setX(1,size[1]*sinGamma);
        latticeVectors[2].setX(0,size[2]*cosBeta);
        latticeVectors[2].setX(1,size[2]*(cosAlpha-cosBeta*cosGamma)/sinGamma);
        latticeVectors[2].setX(2,size[2]*Math.sqrt(1.0-cosAlpha*cosAlpha-cosBeta*cosBeta-cosGamma*cosGamma+2*cosAlpha*cosBeta*cosGamma)/sinGamma);
    }
    
    public void setAngleAlpha(double t) {
        if (t == angle[0]) {
            return;
        }
        setAngles(new double[]{t, angle[1], angle[2]});
    }
    public double getAngleAlpha() {return angle[0];}
    
    public void setBeta(double t) {
        if (t == angle[1]) {
            return;
        }
        setAngles(new double[]{angle[0], t, angle[2]});
    }
    public double getAngleBeta() {return angle[1];}
    
    public void setAngleGamma(double t) {
        if (t == angle[2]) {
            return;
        }
        setAngles(new double[]{angle[0], angle[1], t});
    }
    public double getAngleGamma() {return angle[2];}
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveTriclinic(space, size[0], size[1], size[2], angle[0], angle[1], angle[2]);
    }
        
    public void scaleSize(double scale) {
        setSize(new double[]{size[0]*scale, size[1]*scale, size[2]*scale});
    }        
    
    public int[] latticeIndex(IVector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(IVector q, int[] dimensions) {
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
    
    public String toString() {return "Triclinic";}

}
