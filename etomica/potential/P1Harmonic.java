package etomica.potential;

import etomica.Atom;
import etomica.AtomSet;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Vector;

/**
 * Potential in which attaches a harmonic spring between each affected atom and
 * the nearest boundary in each direction.
 *
 * This class has not been used or checked for correctness.
 *
 * @author David Kofke
 */
 
public class P1Harmonic extends Potential1 implements PotentialSoft {
    
    private final int D;
    private double w = 100.0;
    private final Vector force;
    
    public P1Harmonic(Space space) {
        super(space);
        D = space.D();
        force = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Harmonic potential at the phase boundaries");
        return info;
    }

    public void setSpringConstant(double springConstant) {w = springConstant;}
    public double getSpringConstant() {return w;}
    /**
     * Not implemented correctly.  
     * Returns dimensionless for spring constant.  Should be energy/length^2.
     */
    public etomica.units.Dimension getSpringConstantDimension() {
        return etomica.units.Dimension.NULL;
    }

    public double energy(AtomSet a) {
        Vector r = ((Atom)a).coord.position();
        Vector d = ((Atom)a).node.parentPhase().boundary().dimensions();
        double aSum = 0.0;
        for(int i=0; i<D; i++) {
            double x = r.x(i);
            double dx = d.x(i);
            if(x > 0.5*dx) x = dx - x;
            aSum += x*x;
        }
        return 0.5*w*aSum;
    }
    
    //XXX consider whether 1-body potentials should contribute to virial
    public double virial(AtomSet a) {
        return 0.0;
    }

    public Vector gradient(AtomSet a){
        Vector r = ((Atom)a).coord.position();
        Vector d = ((Atom)a).node.parentPhase().boundary().dimensions();
        force.E(0.0);
        for(int i=0; i<D; i++) {
            double x = r.x(i);
            double dx = d.x(i);
            if(x > 0.5*dx) {
                x = dx - x;
                force.setX(i, -w*x);
            }
            else force.setX(i, w*x);
        }
        return force;
    }//end of gradient
        
}//end of P1Harmonic
   
