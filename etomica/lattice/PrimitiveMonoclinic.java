package etomica.lattice;
import etomica.Space;
import etomica.math.geometry.Polytope;

/**
 * Primitive group for a monoclinic system.  One primitive-vector
 * angle not normal, and vectors not necessarily of equal length.
 * a != b != c; alpha = gamma = 90deg, beta >= 90deg
 */
public class PrimitiveMonoclinic extends Primitive implements Primitive3D {
    
    private boolean isReciprocal = false;
    private double a = 1.0, b = 1.0, c = 1.0;
    private double beta = 0.5*Math.PI, sinBeta = Math.sin(beta), cosBeta = Math.cos(beta);
    
    public PrimitiveMonoclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle);
    }
    public PrimitiveMonoclinic(Space space, double a, double b, double c, double beta) {
        super(space);//also makes reciprocal
        setA(a);//also sets reciprocal via update
        setB(b);
        setC(c);
        setBeta(beta);
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveMonoclinic(Space space, Primitive direct) {
        super(space, direct);
        isReciprocal = true;
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveMonoclinic(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveMonoclinic)reciprocal()).setA(2.0*Math.PI/(a*sinBeta));
        ((PrimitiveMonoclinic)reciprocal()).setB(2.0*Math.PI/b);
        ((PrimitiveMonoclinic)reciprocal()).setC(2.0*Math.PI/(c*sinBeta));
        ((PrimitiveMonoclinic)reciprocal()).setBeta(beta);
    }
    
        //direct lattice (ix = 0, iz = 2)
        // v[0] = (1,0,0); v[1] = (0,1,0); v[2] = (c,0,s)  (times a, b, c)
         
        //reciprocal lattice (ix = 2, iz = 0)
        // v[0] = (s,0,-c); v[1] = (0,1,0); v[2] = (0,0,1);  (times a, b, c)
    public void setA(double a) {
        if(immutable || a <= 0.0) return;
        this.a = a;
        size[0] = a;
        if(!isReciprocal) latticeVectors[0].setX(0,a);
        else {
            latticeVectors[0].setX(0,+a*sinBeta);
            latticeVectors[0].setX(2,-a*cosBeta);
        }
        update();
    }
    public double getA() {return a;}
    
    public void setB(double b) {
        if(immutable || b <= 0.0) return;
        this.b = b;;
        size[1] = b;
        latticeVectors[1].setX(1,b);
        update();
    }
    public double getB() {return b;}
        
    public void setC(double c) {
        if(immutable || c <= 0.0) return;
        this.c = c;
        size[2] = c;
        if(isReciprocal) latticeVectors[2].setX(2,c);
        else {
            latticeVectors[2].setX(0,c*cosBeta);
            latticeVectors[2].setX(2,c*sinBeta);
        }
        update();
    }
    public double getC() {return c;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {
        if(immutable) return;
        if(t < rightAngle) t = rightAngle;
        if(t > Math.PI) t = Math.PI;
        beta = t;
        angle[1] = beta;
        cosBeta = Math.cos(beta);
        sinBeta = Math.sin(beta);
        if(isReciprocal) setA(a); 
        else setC(c);//setA or setC calls update
    }
    public double getBeta() {return beta;}
    
    public void setGamma(double t) {}
    public double getGamma() {return rightAngle;}
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return true;}
    public boolean isEditableC() {return true;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return true;}
    public boolean isEditableGamma() {return false;}


    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveMonoclinic(space, a, b, c, beta);
    }
        
    public void scaleSize(double scale) {
        setA(a*scale);
        setB(b*scale);
        setC(c*scale);
    }        
    
    public int[] latticeIndex(Space.Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
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

}//end of PrimitiveOrthorhombic
    
