/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
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
    protected final Vector dr;
    protected Boundary boundary;
    protected double temperature;
    protected double mass;
    protected double fac, stepFactor = 0.5;
    protected int nPoints = 2;
    
    public P2DiscreteFeynmanHibbs(Space space, Potential2Spherical p2Classical) {
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
        calcFacs();
    }
    
    protected void calcFacs() {
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = stepFactor/Math.sqrt(6*mass/2*temperature/(hbar*hbar));
    }
    
    public void setTemperature(double temperature) {
        this.temperature = temperature;
        calcFacs();
    }
    
    /**
     * Sets the mass; we assume the reduced mass is m/2 (correct for particles
     * with identical mass).
     */
    public void setMass(double m) {
        mass = m;
        calcFacs();
    }

    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }

    public double u(double r2) {
        double r = Math.sqrt(r2);
        if (r < nPoints*fac) {
            return Double.POSITIVE_INFINITY;
        }
        double ueff = p2Classy.u(r2)*r;
        double pnorm = r;
        for (int i=1; i<=nPoints; i++) {
            double pi = Math.exp(-(i*i*stepFactor*stepFactor));
            double ri = r + (i*fac);
            pnorm += pi*ri;
            ueff += p2Classy.u(ri*ri)*pi*ri;
            ri = r - (i*fac);
            pnorm += pi*ri;
            ueff += p2Classy.u(ri*ri)*pi*ri;
        }
        return ueff/pnorm;
    }

    public double getRange() {
        return p2Classy.getRange();
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(1);
        final P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2DiscreteFeynmanHibbs p2dfh = new P2DiscreteFeynmanHibbs(space, p2);
        P2EffectiveFeynmanHibbs p2efh = new P2EffectiveFeynmanHibbs(space, p2);
        double heMass = 4.002602;
        p2dfh.setMass(heMass);
        p2dfh.setTemperature(temperature);
        p2dfh.setStepFactor(0.5);
        p2dfh.setNPoints(2);
        p2efh.setMass(heMass);
        p2efh.setTemperature(temperature);
        for (int i=75;i<5000; i+=50) {
            double r = i/1000.0;
            double uc = p2.u(r*r);
            double ud = p2dfh.u(r*r);
            double ue = p2efh.u(r*r);
        }
    }
}
