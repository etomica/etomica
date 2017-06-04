/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;

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
     * The derivative of the pair energy, times the separation r: r du/dr.
     */
    public abstract double du(double r2);
        
    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public abstract double d2u(double r2);
        
    /**
     * Integral of the potential, used to evaluate corrections for potential truncation.
     * Specifically, this is the integral from rC (the argument) to infinity of
     * u(r) A r^(D-1), where D is the spatial dimension, and A is the area of a unit
     * sphere in D dimensions.  Normally, the long-range potential correction would be obtained
     * by multiplying this quantity by the pair density nPairs/V, where nPairs is the number of pairs of atoms
     * affected by this potential, and V is the volume they occupy.
     */
    public abstract double uInt(double rC);
    
    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        return du(dr.squared());
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return d2u(r2) + du(r2);
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector[] gradient(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < 1.e-10) {
            gradient[0].E(0);
            gradient[1].E(0);
            return gradient;
        }
        gradient[1].Ea1Tv1(du(r2)/r2,dr);
        gradient[0].Ea1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        gradient(atoms);
        pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }
    
    /**
     * Same as uInt.
     */
    public double integral(double rC) {
        return uInt(rC);
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
    protected IBoundary boundary;
    protected final Vector dr;
    
}//end of Potential2SoftSpherical
