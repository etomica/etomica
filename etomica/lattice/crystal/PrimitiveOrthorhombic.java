package etomica.lattice.crystal;
import etomica.Space;
import etomica.lattice.Primitive;
import etomica.lattice.Primitive3D;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Polytope;

/**
 * Primitive group for an orthorhombic system.  All primitive
 * vectors orthogonal but not necessarily of equal length.
 * a != b != c; alpha = beta = gamma = 90deg.
 */
public class PrimitiveOrthorhombic extends Primitive implements Primitive3D {
    
    private double[] sizeCopy;
//    private double a, b, c;
    
    public PrimitiveOrthorhombic(Space space) {
        this(space, 1.0, 1.0, 1.0);
    }
    public PrimitiveOrthorhombic(Space space, double a, double b, double c) {
        super(space); //also makes reciprocal
        //set up orthogonal vectors of unit size
        setA(a);  
        setB(b);
        setC(c); //also sets reciprocal via update
        sizeCopy = new double[space.D];
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveOrthorhombic(Space space, Primitive direct) {
        super(space, direct);
        for(int i=0; i<D; i++) if(size[i] == 0.0) size[i] = 1.0;
//        sizeCopy = new double[space.D()];
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        for(int i=0; i<D; i++) if(size[i] == 0.0) size[i] = 1.0;
        return new PrimitiveOrthorhombic(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveOrthorhombic)reciprocal()).setA(2.0*Math.PI/size[0]);
        ((PrimitiveOrthorhombic)reciprocal()).setB(2.0*Math.PI/size[1]);
        ((PrimitiveOrthorhombic)reciprocal()).setC(2.0*Math.PI/size[2]);
    }
    
    public void setA(double a) {
        if(immutable || a <= 0.0) return;
        size[0] = a;
        latticeVectors[0].setX(0,a);
        update();
    }
    public double getA() {return size[0];}
    
    public void setB(double b) {
        if(immutable || b <= 0.0) return;
        size[1] = b;
        latticeVectors[1].setX(1,b);
        update();
    }
    public double getB() {return size[1];}
        
    public void setC(double c) {
        if(immutable || c <= 0.0) return;
        size[2] = c;
        latticeVectors[2].setX(2,c);
        update();
    }
    public double getC() {return size[2];}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {}
    public double getBeta() {return rightAngle;}
    
    public void setGamma(double t) {}
    public double getGamma() {return rightAngle;}
    

    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return true;}
    public boolean isEditableC() {return true;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return false;}
    public boolean isEditableGamma() {return false;}
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveOrthorhombic(space, size[0], size[1], size[2]);
    }
    
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSize(double size) {
        if(immutable) return;
        for(int i=0; i<D; i++) {
            this.size[i] = size;
            latticeVectors[i].setX(i,size);
        }
        update();
    }
    
    /**
     * Returns a copy of the array of primitive-vector sizes.
     */
     //used by IteratorFactoryCell
    public double[] getSize() {
        for(int i=0; i<D; i++) sizeCopy[i] = size[i];
        return sizeCopy;
    }
    
    public void scaleSize(double scale) {
        if(immutable) return;
        for(int i=0; i<D; i++) {
            size[i] *= scale;
            latticeVectors[i].setX(i,size[i]);
        }
        update();
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
    
    /**
     * Returns a new Cuboid with edges of length given by the current
     * values of the primitive vectors.
     */
    public Polytope unitCell() {
        return new Cuboid(size[0], size[1], size[2]);
    }
    
    public String toString() {return "Orthorhombic";}

}//end of PrimitiveOrthorhombic
    
