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
import etomica.units.Undefined;

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
    private final Vector x0;
    
    public P1Harmonic(Space space) {
        super(space);
        D = space.D();
        force = space.makeVector();
        x0 = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Harmonic potential at the phase boundaries");
        return info;
    }

    public void setSpringConstant(double springConstant) {
        w = springConstant;
        System.out.println("spring constant " + w);
    }
    
    public double getSpringConstant() {
        return w;
    }
    
    public void setX0(Vector x0) {
        this.x0.E(x0);
        System.out.println(x0);
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
        return 0.5*w*((AtomLeaf)a).coord.position().Mv1Squared(x0);
//        r.Mv1Squared(x0);
//        Vector d = ((Atom)a).node.parentPhase().boundary().dimensions();
//        double aSum = 0.0;
//        for(int i=0; i<D; i++) {
//            double x = r.x(i) - x0;
//            double dx = d.x(i);
//            if(x > 0.5*dx) x = dx - x;
//            aSum += x*x;
//        }
//        return 0.5*w*aSum;
    }
    
    //XXX consider whether 1-body potentials should contribute to virial
    public double virial(AtomSet a) {
        return 0.0;
    }

    public Vector gradient(AtomSet a){
        Vector r = ((AtomLeaf)a).coord.position();
//        Vector d = ((Atom)a).node.parentPhase().boundary().dimensions();
        force.Ev1Mv2(r,x0);
        force.TE(w);
//        for(int i=0; i<D; i++) {
//            double x = r.x(i) - x0;
//            double dx = d.x(i);
//            if(x > 0.5*dx) {
//                x = dx - x;
//                force.setX(i, -w*x);
//            }
//            else force.setX(i, w*x);
//        }
            
        return force;
    }//end of gradient
        
}//end of P1Harmonic
   
