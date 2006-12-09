package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.space.Space;
import etomica.space.Vector;
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
    private final Vector[] force;
    private final Vector x0;
    
    public P1Harmonic(Space space) {
        super(space);
        force = new Vector[]{space.makeVector()};
        x0 = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Harmonic potential at the phase boundaries");
        return info;
    }

    public void setSpringConstant(double springConstant) {
        w = springConstant;
    }
    
    public double getSpringConstant() {
        return w;
    }
    
    public void setX0(Vector x0) {
        this.x0.E(x0);
    }
    
    public Vector getX0() {
        return (Vector)x0.clone();
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

    public double energy(AtomSet a) {
        return 0.5*w*((AtomLeaf)a).getCoord().position().Mv1Squared(x0);
    }
    
    //XXX consider whether 1-body potentials should contribute to virial
    public double virial(AtomSet a) {
        return 0.0;
    }

    public Vector[] gradient(AtomSet a){
        Vector r = ((AtomLeaf)a).getCoord().position();
        force[0].Ev1Mv2(r,x0);
        force[0].TE(w);
            
        return force;
    }
        
}
   
