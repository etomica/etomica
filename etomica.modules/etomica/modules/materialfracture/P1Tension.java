package etomica.modules.materialfracture;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Null;

public class P1Tension implements PotentialSoft {
    
    protected final ISpace space;
    protected double w;
    protected final IVector[] force;
    protected IBox box;
    
    public P1Tension(ISpace space) {
        this.space = space;
        force = new IVector[1];
        force[0] = space.makeVector();
        setSpringConstant(100.0);
    }
    
    public int nBody() {
        return 1;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public void setSpringConstant(double springConstant) {w = springConstant;}
    public double getSpringConstant() {return w;}
    
    /**
     * Not implemented correctly.  
     * Returns dimensionless for spring constant.  Should be energy/length^2.
     */
    public etomica.units.Dimension getSpringConstantDimension() {
        return Null.DIMENSION;
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }

    public double energy(IAtomSet a) {
        IVector r = ((IAtomPositioned)a.getAtom(0)).getPosition();
        double aSum = 0.0;
        double x = r.x(0);
        aSum += x*x;
        return 0.5*w*aSum;
    }

    public IVector[] gradient(IAtomSet a, Tensor t){
        return gradient(a);
    }
    
    public IVector[] gradient(IAtomSet a) {
        IVector r = ((IAtomPositioned)a.getAtom(0)).getPosition();
        force[0].setX(1, 0);
        double x = r.x(0);
        force[0].setX(0,-w*x);
        return force;
    }
    
    public double virial (IAtomSet a) {
        return 0;
    }
}
   
