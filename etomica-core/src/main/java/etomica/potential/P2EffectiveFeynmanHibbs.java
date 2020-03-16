/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Hartree;
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
    protected final Vector dr;
    protected Boundary boundary;
    protected double temperature;
    protected double mass;
    protected double fac;
    
    public P2EffectiveFeynmanHibbs(Space space, Potential2SoftSpherical p2Classical) {
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
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }

    public double u(double r2) {
        double uc = p2Classy.u(r2);
        if (Double.isInfinite(uc)) { return uc; }
        double duc = p2Classy.du(r2); //multiplied by r
        double d2uc = p2Classy.d2u(r2);  //multiplied by r*r
        if (d2uc == Double.POSITIVE_INFINITY) return d2uc;
        double u = uc + fac*(d2uc + 2*duc)/r2;
        // if the classical potential is repulsive, the semiclassical potential
        // should be more repulsive.  In nearly all cases, it is, but this often
        // fails for very short separations.  Just enforce it here.
        if (uc > 0 && u < uc) return uc;
        return u;
    }
    
    public double dudkT(double r2) {
        double uc = p2Classy.u(r2);
        if (Double.isInfinite(uc)) { return 0; }
        double duc = p2Classy.du(r2);
        double d2uc = p2Classy.d2u(r2);
        if (Double.isInfinite(fac*(d2uc + 2*duc)/r2*(-2/temperature/temperature))) {
            throw new RuntimeException("fac*(d2uc + 2*duc)/r2 is infinite");
        }

        return fac * (d2uc + 2 * duc) / r2 / (-temperature);  //fac includes temperature in the denominator
    }

    public void setBox(Box box) {
        p2Classy.setBox(box);
        boundary = box.getBoundary();
    }
    
    public double getRange() {
        return p2Classy.getRange();
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(20);
        final P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2EffectiveFeynmanHibbs p2fh = new P2EffectiveFeynmanHibbs(space, p2);
        double heMass = 4.002602;
        p2fh.setMass(heMass);
        p2fh.setTemperature(temperature);
        for (int i=22;i<2000; i++) {
            double r = i/100.0;
            double u = Hartree.UNIT.fromSim(p2fh.u(r*r));
            System.out.println(r+" "+u);
        }
        

    }
}
