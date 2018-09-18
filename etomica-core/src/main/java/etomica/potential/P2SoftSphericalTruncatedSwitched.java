/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  The energy
 * is switched from fully-on to 0 over a short range (i.e. from 0.95*rC to rC).
 */
public class P2SoftSphericalTruncatedSwitched extends Potential2 implements PotentialSoft {
    
    public P2SoftSphericalTruncatedSwitched(Space _space, Potential2SoftSpherical potential,
                                            double truncationRadius) {
        super(_space);
        this.potential = potential;
        setTruncationRadius(truncationRadius);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        setSwitchFac(0.95);
    }

    /**
     * Returns the wrapped potential.
     */
    public PotentialSoft getWrappedPotential() {
        return potential;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
        r2Switch = r2Cutoff*switchFac*switchFac;
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {return rCutoff;}
    
    public double getSwitchFac() {
        return switchFac;
    }

    public void setSwitchFac(double newSwitchFac) {
        switchFac = newSwitchFac;
        r2Switch = r2Cutoff*switchFac*switchFac;
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }

    public Vector[] gradient(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < r2Cutoff) {
            Vector[] unswitchedGradient = potential.gradient(atoms);
            gradient[0].E(unswitchedGradient[0]);
            gradient[1].E(unswitchedGradient[1]);
            if (r2 > r2Switch) {
                double r = Math.sqrt(r2);
                double fac = getF(r);
                double u = potential.u(r2);
                gradient[0].TE(fac);
                gradient[1].TE(fac);
                // (df/dr)/r
                fac = getdFdr(r)/r;
                gradient[0].PEa1Tv1(-fac*u, dr);
                gradient[1].PEa1Tv1(+fac*u, dr);
            }
            return gradient;
        }
        gradient[0].E(0);
        gradient[1].E(0);
        return gradient;
    }
    
    protected double getF(double r) {
        switch (taperOrder) {
            case 1:
                return (rCutoff-r)/(rCutoff*(1-switchFac));
            case 2:
                return (r2Cutoff-2*rCutoff*r+r*r)/(r2Cutoff*(1-switchFac)*(1-switchFac))+1e-7;
            case 3:
                double rt = switchFac*rCutoff;
                double a = (r-rt)/(rCutoff-rt);
                return (rCutoff-r)/(rCutoff-rt)*(1-a*a) + a*(1-a)*(1-a);
            default:
                throw new RuntimeException("oops");
        }
    }

    protected double getdFdr(double r) {
        switch (taperOrder) {
            case 1:
                return -1.0/(rCutoff*(1-switchFac));
            case 2:
                return -2 * (rCutoff - r) / (r2Cutoff*(1-switchFac)*(1-switchFac));
            case 3:
                double rt = switchFac*rCutoff;
                double a = (r-rt)/(rCutoff-rt);
                double b = rCutoff-rt;
                double c = rCutoff-r;
                return (-(1.0-a*a)/b - 2*c*a/(b*b)) + (1-a)*(1-a)/b - 2*a/b + 2*a*a/b;
            default:
                throw new RuntimeException("oops");
        }
    }
    
    public static void main(String[] args) {
    	Space sp = Space3D.getInstance();
        P2SoftSphericalTruncatedSwitched p = new P2SoftSphericalTruncatedSwitched(sp, new P2LennardJones(sp), 2);
        for (double x = 1.900001; x<1.999999; x+=0.001) {
            System.out.println(x+" "+p.getF(x)+" "+p.getdFdr(x));
        }
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
    
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (dr.squared() > r2Cutoff) {
            return 0;
        }
        double u = potential.u(r2);
        if (dr.squared() > r2Switch) {
            u *= getF(Math.sqrt(r2));
        }
        return u;
    }
    
    public double virial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        if (dr.squared() < r2Cutoff) {
            return potential.virial(atoms);
        }
        return 0;
    }
    
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {return Length.DIMENSION;}
    
    public void setBox(Box newBox) {
        potential.setBox(newBox);
        boundary = newBox.getBoundary();
    }
    
    private static final long serialVersionUID = 1L;
    protected double rCutoff, r2Cutoff;
    protected final Potential2SoftSpherical potential;
    protected final Vector dr;
    protected Boundary boundary;
    protected final Vector[] gradient;
    protected int taperOrder = 3;
    protected double switchFac, r2Switch;
}
