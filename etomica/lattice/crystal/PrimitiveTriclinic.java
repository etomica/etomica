package etomica.lattice.crystal;
import etomica.Space;
import etomica.Space3D;
import etomica.lattice.Primitive;
import etomica.math.geometry.Polytope;
import etomica.space.Vector;

/**
 * Primitive group for a triclinic system.  No restrictions on
 * primitive-vector angles or lengths.
 * a != b != c; alpha != gamma != beta
 */
public class PrimitiveTriclinic extends Primitive implements Primitive3D {
    
    private double a = 1.0, b = 1.0, c = 1.0;
    private double alpha = 0.5*Math.PI, sinAlpha = Math.sin(alpha), cosAlpha = Math.cos(alpha);
    private double beta = 0.5*Math.PI, sinBeta = Math.sin(beta), cosBeta = Math.cos(beta);
    private double gamma = 0.5*Math.PI, sinGamma = Math.sin(gamma), cosGamma = Math.cos(gamma);
    private boolean isReciprocal = false;
    
    public PrimitiveTriclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle, rightAngle, rightAngle);
    }
    public PrimitiveTriclinic(Space space, double a, double b, double c, 
                                              double alpha, double beta, double gamma) {
        super(space);
        setA(a);
        setB(b);
        setC(c);
        setAlpha(alpha);
        setBeta(beta);
        setGamma(gamma);
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveTriclinic(Space space, Primitive direct) {
        super(space, direct);
        isReciprocal = true;
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveTriclinic(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        PrimitiveTriclinic recip = (PrimitiveTriclinic)reciprocal();
        Vector aStar = (Vector)recip.latticeVectors[0];
        Vector bStar = (Vector)recip.latticeVectors[1];
        Vector cStar = (Vector)recip.latticeVectors[2];
        Vector aVec = (Vector)latticeVectors[0];
        Vector bVec = (Vector)latticeVectors[1];
        Vector cVec = (Vector)latticeVectors[2];
        aStar.E(bVec);
        aStar.XE(cVec);
        double factor = 2.0*Math.PI/aVec.dot(aStar); // a . (b X c)
        aStar.TE(factor);
        bStar.E(cVec);
        bStar.XE(aVec);
        bStar.TE(factor);
        cStar.E(aVec);
        cStar.XE(bVec);
        cStar.TE(factor);
    }
    
    public void setA(double a) {
        if(immutable || a <= 0.0) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        this.a = a;
        size[0] = a;
        latticeVectors[0].setX(0,a);
        update();
    }
    public double getA() {return a;}
    
    public void setB(double b) {
        if(immutable || b <= 0.0) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        this.b = b;
        size[1] = b;
        latticeVectors[1].setX(0,b*cosGamma);
        latticeVectors[1].setX(1,b*sinGamma);
        update();
    }
    public double getB() {return b;}
        
    public void setC(double c) {
        if(immutable || c <= 0.0) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        this.c = c;
        size[2] = c;
        latticeVectors[2].setX(0,c*cosBeta);
        latticeVectors[2].setX(1,c*(cosAlpha-cosBeta*cosGamma)/sinGamma);
        latticeVectors[2].setX(2,c*Math.sqrt(1.0-cosAlpha*cosAlpha-cosBeta*cosBeta-cosGamma*cosGamma+2*cosAlpha*cosBeta*cosGamma)/sinGamma);
        update();
    }
    public double getC() {return c;}
    
    private double bounds(double t) {
        if(t < 0.0) return 0.0;
        else if(t > Math.PI) return Math.PI;
        else return t;
    }
    
    public void setAlpha(double t) {
        if(immutable) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        t = bounds(t);
        alpha = t;
        angle[0] = alpha;
        cosAlpha = Math.cos(alpha);
        sinAlpha = Math.sin(alpha);
        setC(c);
    }
    public double getAlpha() {return alpha;}
    
    public void setBeta(double t) {
        if(immutable) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        t = bounds(t);
        beta = t;
        angle[1] = beta;
        cosBeta = Math.cos(beta);
        sinBeta = Math.sin(beta);
        setC(c);
    }
    public double getBeta() {return beta;}
    
    public void setGamma(double t) {
        if(immutable) return;
        if(isReciprocal) throw new RuntimeException("Error: PrimitiveTriclinic reciprocal vectors cannot be edited directly");
        t = bounds(t);
        gamma = t;
        angle[2] = gamma;
        cosGamma = Math.cos(gamma);
        sinGamma = Math.sin(gamma);
        setB(b);
        setC(c);
    }
    public double getGamma() {return gamma;}
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return true;}
    public boolean isEditableC() {return true;}
    public boolean isEditableAlpha() {return true;}
    public boolean isEditableBeta() {return true;}
    public boolean isEditableGamma() {return true;}

    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveTriclinic(space, a, b, c, alpha, beta, gamma);
    }
        
    public void scaleSize(double scale) {
        setA(a*scale);
        setB(b*scale);
        setC(c*scale);
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
        throw new RuntimeException("method PrimitiveOrthorhombic.unitCell not yet implemented");
    }
    
    public String toString() {return "Triclinic";}

}//end of PrimitiveOrthorhombic
    
