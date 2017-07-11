/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.box.Box;
import etomica.models.water.P2WaterSPCSoft;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeOrientedDynamic;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  The energy
 * is switched from fully-on to 0 over a short range (i.e. from 0.95*rC to rC).
 */
public class P2MoleculeSoftTruncatedSwitched extends PotentialMolecular implements IPotentialMolecularTorque {
    
    public P2MoleculeSoftTruncatedSwitched(IPotentialMolecularTorque potential, double truncationRadius, Space _space) {
        super(2, _space);
        this.potential = potential;
        setTruncationRadius(truncationRadius);
        gradientAndTorque = new Vector[2][2];
        for (int i=0; i<2; i++) {
            gradientAndTorque[0][i] = space.makeVector();
            gradientAndTorque[1][i] = space.makeVector();
        }
        dr = space.makeVector();
        setSwitchFac(0.95);
    }
    
    /**
     * Returns the wrapped potential.
     */
    public IPotentialMolecularTorque getWrappedPotential() {
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

    public double getSwitchFac() {
        return switchFac;
    }

    public void setSwitchFac(double newSwitchFac) {
        switchFac = newSwitchFac;
        r2Switch = r2Cutoff*switchFac*switchFac;
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {return rCutoff;}
    
    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }

    public Vector[][] gradientAndTorque(IMoleculeList atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getMolecule(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getMolecule(0)).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < r2Cutoff) {
            Vector[][] gradientAndTorqueUnswitched = potential.gradientAndTorque(atoms);
            gradientAndTorque[0][0].E(gradientAndTorqueUnswitched[0][0]);
            gradientAndTorque[0][1].E(gradientAndTorqueUnswitched[0][0]);
            gradientAndTorque[1][0].E(gradientAndTorqueUnswitched[1][0]);
            gradientAndTorque[1][1].E(gradientAndTorqueUnswitched[1][0]);
            if (r2 > r2Switch) {
                double r = Math.sqrt(r2);
                double fac = getF(r);
                // G = f G + u df/dr
                // tau = f tau
                gradientAndTorque[0][0].TE(fac);
                gradientAndTorque[0][1].TE(fac);
                gradientAndTorque[1][0].TE(fac);
                gradientAndTorque[1][1].TE(fac);
                // (df/dr)/r
                fac = getdFdr(r)/r;
                double u = potential.energy(atoms);
                gradientAndTorque[0][0].PEa1Tv1(-fac*u, dr);
                gradientAndTorque[0][1].PEa1Tv1(+fac*u, dr);
            }
            return gradientAndTorqueUnswitched;
        }
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(0);
        gradientAndTorque[1][1].E(0);
        return gradientAndTorque;
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
        P2MoleculeSoftTruncatedSwitched p095 = new P2MoleculeSoftTruncatedSwitched(new P2WaterSPCSoft(sp), 2, sp);
        p095.setSwitchFac(0.95);
        P2MoleculeSoftTruncatedSwitched p080 = new P2MoleculeSoftTruncatedSwitched(new P2WaterSPCSoft(sp), 2, sp);
        p080.setSwitchFac(0.80);
        P2MoleculeSoftTruncatedSwitched p050 = new P2MoleculeSoftTruncatedSwitched(new P2WaterSPCSoft(sp), 2, sp);
        p050.setSwitchFac(0.50);
        P2MoleculeSoftTruncatedSwitched p010 = new P2MoleculeSoftTruncatedSwitched(new P2WaterSPCSoft(sp), 2, sp);
        p010.setSwitchFac(0.10);
        for (double x = 0.002; x<2; x+=0.002) {
            System.out.println(x+" "+(1.0/x)+" "+p095.getF(x)/x+" "+p080.getF(x)/x+" "+p050.getF(x)/x+" "+p010.getF(x)/x);
        }
    }
    
    public Vector[] gradient(IMoleculeList atoms) {
        return gradientAndTorque(atoms)[0];
    }

    public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
        return gradientAndTorque(atoms)[0];
    }
    
    public double energy(IMoleculeList atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getMolecule(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getMolecule(0)).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (dr.squared() > r2Cutoff) {
            return 0;
        }
        double u = potential.energy(atoms);
        if (r2 > r2Switch) {
            u *= getF(Math.sqrt(r2));
        }
        return u;
    }
    
    public double virial(IMoleculeList atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getMolecule(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getMolecule(0)).getPosition());
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
    protected final IPotentialMolecularTorque potential;
    protected final Vector dr;
    protected Boundary boundary;
    protected final Vector[][] gradientAndTorque;
    protected int taperOrder = 3;
    protected double switchFac, r2Switch;
}
