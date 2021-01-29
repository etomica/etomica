/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public abstract class Potential2SoftSpherical extends Potential2 implements Potential2Soft {

    public Potential2SoftSpherical(Space space) {
        super(space);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
    }

    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        return du(dr.squared());
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return d2u(r2) + du(r2);
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector[] gradient(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double du = du(r2);
        if (du == 0) {
            gradient[0].E(0.0);
            gradient[1].E(0.0);
            return gradient;
        }
        gradient[1].Ea1Tv1(du/r2,dr);
        gradient[0].Ea1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        gradient(atoms);
        if (!gradient[0].isZero()) {
            pressureTensor.PEv1v2(gradient[0],dr);
        }
        return gradient;
    }

    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    protected final Vector[] gradient;
    protected Boundary boundary;
    protected final Vector dr;
    
}//end of Potential2SoftSpherical
