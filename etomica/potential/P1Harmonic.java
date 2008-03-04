package etomica.potential;

import etomica.EtomicaInfo;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Potential in which attaches a harmonic spring between each affected atom and
 * the nearest boundary in each direction.
 *
 * This class has not been used or checked for correctness.
 *
 * @author David Kofke
 */
 
public class P1Harmonic extends Potential1 implements PotentialSoft {
    
    private static final long serialVersionUID = 1L;
    private double w = 100.0;
    private final IVector[] force;
    private final IVector x0;
    
    public P1Harmonic(Space space) {
        super(space);
        force = new IVector[]{space.makeVector()};
        x0 = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Harmonic potential at the box boundaries");
        return info;
    }

    public void setSpringConstant(double springConstant) {
        w = springConstant;
    }
    
    public double getSpringConstant() {
        return w;
    }
    
    public void setX0(IVector x0) {
        this.x0.E(x0);
    }
    
    public IVector getX0() {
        return x0;
    }
    
    public Dimension getX0Dimension() {
        return Length.DIMENSION;
    }
    /**
     * Not implemented correctly.  
     * Should be energy/length^2.
     */
    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-2});
    }

    public double energy(IAtomSet a) {
        return 0.5*w*((IAtomPositioned)a.getAtom(0)).getPosition().Mv1Squared(x0);
    }
    
    public double virial(IAtomSet a) {
        return 0.0;
    }

    public IVector[] gradient(IAtomSet a){
        IVector r = ((IAtomPositioned)a.getAtom(0)).getPosition();
        force[0].Ev1Mv2(r,x0);
        force[0].TE(w);
            
        return force;
    }
        
    public IVector[] gradient(IAtomSet a, Tensor pressureTensor){
        return gradient(a);
    }
}
   
