/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncated extends Potential2SoftSpherical
               implements PotentialTruncated {

    protected final Potential2SoftSpherical potential;
    protected double rCutoff, r2Cutoff;
    protected boolean makeLrc = true;

    public P2SoftSphericalTruncated(Space _space, Potential2SoftSpherical potential, double truncationRadius) {
        super(_space);
        this.potential = potential;
        setTruncationRadius(truncationRadius);
    }

    public P2SoftSphericalTruncated(Potential2SoftSpherical potential, double truncationRadius) {
        this(potential.space, potential, truncationRadius);
    }

    /**
     * Returns the wrapped potential.
     */
    public Potential2SoftSpherical getWrappedPotential() {
        return potential;
    }

    public void setBox(Box box) {
        potential.setBox(box);
        super.setBox(box);
    }

    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double u(double r2) {
        return (r2 < r2Cutoff) ? potential.u(r2) : 0.0;
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double du(double r2) {
        return (r2 < r2Cutoff) ? potential.du(r2) : 0.0;
    }

    /**
     * Returns the 2nd derivative (r^2 d^2u/dr^2) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double d2u(double r2) {
        return (r2 < r2Cutoff) ? potential.d2u(r2) : 0.0;
    }

    /**
     * Returns the value of uInt for the wrapped potential.
     */
    public double uInt(double rC) {
        return potential.uInt(rC);
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {
        return rCutoff;
    }
    
    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
    }
    
    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }
    
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {return Length.DIMENSION;}
    
    /**
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(AtomType[] types) {
        if (!makeLrc) return null;
        return new P0Lrc(space, potential, this, types);
    }
    
    public boolean getMakeLrc() {
        return makeLrc;
    }

    public void setMakeLrc(boolean newMakeLrc) {
        makeLrc = newMakeLrc;
    }

    /**
     * Inner class that implements the long-range correction for this truncation scheme.
     */
    public static class P0Lrc extends Potential0Lrc implements Potential2Soft {

        private final double A;
        private final int D;
        private Potential2Soft potential;

        public P0Lrc(Space space, Potential2Soft truncatedPotential,
                     Potential2Soft potential, AtomType[] types) {
            super(space, types, truncatedPotential);
            this.potential = potential;
            A = space.sphereArea(1.0);  //multiplier for differential surface element
            D = space.D();              //spatial dimension
        }

        public double energy(IAtomList atoms) {
            return uCorrection(nPairs()/box.getBoundary().volume());
        }

        public double virial(IAtomList atoms) {
            return duCorrection(nPairs()/box.getBoundary().volume());
        }

        public double hyperVirial(IAtomList pair) {
            return d2uCorrection(nPairs()/box.getBoundary().volume()) + duCorrection(nPairs()/box.getBoundary().volume());
        }

        public Vector[] gradient(IAtomList atoms) {
            return null;
        }

        public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
            double virial = virial(atoms) / pressureTensor.D();
            for (int i=0; i<pressureTensor.D(); i++) {
                pressureTensor.setComponent(i,i,pressureTensor.component(i,i)-virial);
            }
            // we'd like to throw an exception and return the tensor, but we can't.  return null
            // instead.  it should work about as well as throwing an exception.
            return null;
        }

        /**
         * Uses result from integration-by-parts to evaluate integral of
         * r2 d2u/dr2 using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
        public double d2uCorrection(double pairDensity) {
            double rCutoff = potential.getRange();
            double integral = ((Potential2Soft)truncatedPotential).integral(rCutoff);
            //need potential to be spherical to apply here
            integral = -A*space.powerD(rCutoff)*((Potential2Soft)truncatedPotential).u(rCutoff*rCutoff) - D*integral;
            integral = -A*space.powerD(rCutoff)*((Potential2SoftSpherical)truncatedPotential).du(rCutoff*rCutoff) - (D+1)*integral;
            return pairDensity*integral;
        }

        /**
         * Long-range correction to the energy.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
        public double uCorrection(double pairDensity) {
            double rCutoff = potential.getRange();
            double integral = ((Potential2Soft)truncatedPotential).integral(rCutoff);
            return pairDensity*integral;
        }

        /**
         * Uses result from integration-by-parts to evaluate integral of
         * r du/dr using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
        public double duCorrection(double pairDensity) {
            double rCutoff = potential.getRange();
            double integral = ((Potential2Soft)truncatedPotential).integral(rCutoff);
            //need potential to be spherical to apply here
            integral = -A*space.powerD(rCutoff)*((Potential2Soft)truncatedPotential).u(rCutoff*rCutoff) - D*integral;
            return pairDensity*integral;
        }

        public double u(double r2) {
            throw new RuntimeException("nope");
        }

        public double integral(double rC) {
            throw new RuntimeException("nope");
        }

        public double du(double r2) {
            throw new RuntimeException("nope");
        }
    }//end of P0lrc
}
