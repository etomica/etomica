package etomica.virial.simulations;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.potential.P2HePCKLJSFancy;
import etomica.potential.Potential2Spherical;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;

/**
 * Effective 2 body potential approximating the quantum behavior of atomic
 * interactions.
 * 
 * This is based on the same idea as Feynman-Hibbs effective potential
 * (approximating the distribution of atomic centers with a Gaussian), but
 * uses a weighted average of the potential evaluated at discrete points,
 * rather than a Taylor series expansion.  As such, it is less accurate at
 * high temperatures (where the series works well because the probability
 * distribution is narrow) but can work better at low temperature (where
 * the series fails because the distribution is very wide).
 *  
 * @author Andrew Schultz
 */
public class P2DiscreteFeynmanHibbs implements Potential2Spherical {

    protected final Potential2Spherical p2Classy;
    protected final IVectorMutable dr;
    protected IBoundary boundary;
    protected double temperature;
    protected double mass;
    protected double fac, stepFactor = 2.0/3.0;
    protected int nPoints = 2;
    
    public P2DiscreteFeynmanHibbs(ISpace space, Potential2Spherical p2Classical) {
        p2Classy = p2Classical;
        dr = space.makeVector();
    }
    
    public int nBody() {
        return 2;
    }
    
    public void setNPoints(int nPoints) {
        this.nPoints = nPoints;
    }
    
    public void setStepFactor(double stepFactor) {
        this.stepFactor = stepFactor;
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = stepFactor/Math.sqrt(6*mass/2*temperature/(hbar*hbar));
    }
    
    /**
     * Sets the mass; we assume the reduced mass is m/2 (correct for particles
     * with identical mass).
     */
    public void setMass(double m) {
        mass = m;
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = stepFactor/Math.sqrt(6*mass/2*temperature/(hbar*hbar));
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
        double r = Math.sqrt(r2);
        if (r < nPoints*fac) {
            return Double.POSITIVE_INFINITY;
        }
        double ueff = p2Classy.u(r2)*r*r;
        double pnorm = r*r;
        for (int i=1; i<=nPoints; i++) {
            double pi = Math.exp(-(i*i*stepFactor*stepFactor));
            double ri = r + (i*fac);
            pnorm += pi*ri*ri;
            ueff += p2Classy.u(ri*ri)*pi*ri*ri;
            ri = r - (i*fac);
            pnorm += pi*ri*ri;
            ueff += p2Classy.u(ri*ri)*pi*ri*ri;
        }
        return ueff/pnorm;
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
        final P2HePCKLJSFancy p2 = new P2HePCKLJSFancy(space);
        P2DiscreteFeynmanHibbs p2SemiClassical = new P2DiscreteFeynmanHibbs(space, p2);
        double heMass = 4.002602;
        p2SemiClassical.setMass(heMass);
        p2SemiClassical.setTemperature(temperature);
        for (int i=25;i<500; i++) {
            double r = i/100.0;
            double u = p2SemiClassical.u(r*r);
            System.out.println(r+" "+Math.exp(-u/temperature));
        }
    }
}
