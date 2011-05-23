package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;

/**
 * Effective 2 body potential approximating the quantum behavior of atomic
 * interactions.
 *  
 *  R. P. Feynman and A. R. Hibbs, Quantum Mechanics and Path Integrals
 *   McGraw-Hill, New York, 1965 ; R. P. Feynman, Statistical Mechanics:
 *  Guillot B. and Guissani Y., "Quantum effects in simulated water by 
 *   the Feynman­Hibbs approach," J. Chem Phys., 108 (1998) 10162
 *  
 * @author Andrew Schultz
 */
public class P2EffectiveFeynmanHibbs implements Potential2Spherical {

    protected final Potential2SoftSpherical p2Classy;
    protected final IVectorMutable dr;
    protected IBoundary boundary;
    protected double temperature;
    protected double mass;
    protected double fac;
    
    public P2EffectiveFeynmanHibbs(ISpace space, Potential2SoftSpherical p2Classical) {
        p2Classy = p2Classical;
        dr = space.makeVector();
    }
    
    public int nBody() {
        return 2;
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = hbar*hbar/(24*mass/2)/temperature;
    }
    
    /**
     * Sets the mass; we assume the reduced mass is m/2 (correct for particles
     * with identical mass).
     */
    public void setMass(double m) {
        mass = m;
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = hbar*hbar/(24*m/2)/temperature;
    }

    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }

    public double u(double r2) {
        double uc = p2Classy.u(r2);
        double duc = p2Classy.du(r2);
        double d2uc = p2Classy.d2u(r2);
        return uc + fac*(d2uc + 2*duc)/r2;
    }

    public void setBox(IBox box) {
        p2Classy.setBox(box);
        boundary = box.getBoundary();
    }
    
    public double getRange() {
        return p2Classy.getRange();
    }

    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(10);
        final P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2DiscreteFeynmanHibbs p2fh = new P2DiscreteFeynmanHibbs(space, p2);
        double heMass = 4.002602;
        p2fh.setMass(heMass);
        p2fh.setTemperature(temperature);
        for (int i=25;i<1000; i++) {
            double r = i/100.0;
            double u = p2fh.u(r*r);
            System.out.println(r+" "+u);
        }
    }
}
