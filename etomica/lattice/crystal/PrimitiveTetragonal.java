package etomica.lattice.crystal;
import etomica.lattice.Primitive;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Primitive group for a tetragonal system.  All primitive
 * vectors orthogonal and two are of equal length.
 */
public class PrimitiveTetragonal extends Primitive implements Primitive3D {
    
    private double ab = 1.0, c = 1.0;
    
    public PrimitiveTetragonal(Space space) {
        this(space, 1.0, 1.0);
    }
    public PrimitiveTetragonal(Space space, double ab, double c) {
        super(space);//also makes reciprocal
        setAB(ab); //also sets reciprocal via update
        setC(c);
    }
    /**
     * Constructor used by makeReciprocal method.
     */
    private PrimitiveTetragonal(Space space, Primitive direct) {
        super(space, direct);
    }
    
    //called by superclass constructor
    protected Primitive makeReciprocal() {
        return new PrimitiveTetragonal(space, this);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        ((PrimitiveTetragonal)reciprocal()).setAB(2.0*Math.PI/ab);
        ((PrimitiveTetragonal)reciprocal()).setC(2.0*Math.PI/c);
    }
    
    public void setA(double a) {setAB(a);}
    public double getA() {return ab;}
    
    public void setB(double b) {setAB(b);}
    public double getB() {return ab;}
    
    public void setAB(double ab) {
        if(immutable || ab <= 0.0) return;
        this.ab = ab;
        size[0] = size[1] = ab;
        latticeVectors[0].setX(0,ab);
        latticeVectors[1].setX(1,ab);
        update();
    }
    
    public void setC(double c) {
        if(immutable || c <= 0.0) return;
        this.c = c;
        size[2] = c;
        latticeVectors[2].setX(2, c);
        update();
    }
    public double getC() {return c;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {}
    public double getBeta() {return rightAngle;}
    
    public void setGamma(double t) {}
    public double getGamma() {return rightAngle;}
    
    
    public boolean isEditableA() {return true;}
    public boolean isEditableB() {return false;}
    public boolean isEditableC() {return true;}
    public boolean isEditableAlpha() {return false;}
    public boolean isEditableBeta() {return false;}
    public boolean isEditableGamma() {return false;}

    /**
     * Returns a new PrimitiveTetragonal with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveTetragonal(space, ab, c);
    }
    
    public void scaleSize(double scale) {
        setAB(ab*scale);
        setC(c*scale);
    }

    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
   */ }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
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
        return new Cuboid(space, ab, ab, c);
    }
    
    public String toString() {return "Tetragonal";}
}
