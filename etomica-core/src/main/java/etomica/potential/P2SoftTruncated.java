/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Boundary;
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
public class P2SoftTruncated extends Potential2
               implements PotentialTruncated, Potential2Soft {

    private static final long serialVersionUID = 1L;
    protected final Vector dr;
    protected final Potential2Soft wrappedPotential;
    protected final Vector[] gradient;
    protected double rCutoff, r2Cutoff;
    protected Boundary boundary;

    public P2SoftTruncated(Potential2Soft potential, double truncationRadius, Space _space) {
        super(_space);
        this.wrappedPotential = potential;
        setTruncationRadius(truncationRadius);
        dr = space.makeVector();
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
    }
    
    /**
     * Returns the wrapped potential.
     */
    public Potential2Soft getWrappedPotential() {
        return wrappedPotential;
    }

    public void setBox(Box newBox) {
        wrappedPotential.setBox(newBox);
        boundary = newBox.getBoundary();
    }
    
    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return (r2 < r2Cutoff) ? wrappedPotential.energy(atoms) : 0;
    }

    /**
     * Returns the 2nd derivative (r^2 d^2u/dr^2) of the wrapped potential if the separation
     * is less than the cutoff value
     */
    public double virial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return (r2 < r2Cutoff) ? wrappedPotential.virial(atoms) : 0;
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     */
    public Vector[] gradient(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return (r2 < r2Cutoff) ? wrappedPotential.gradient(atoms) : gradient;
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     */
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 > r2Cutoff) return gradient;
        return wrappedPotential.gradient(atoms, pressureTensor);
    }
    
    public double hyperVirial(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return (r2 < r2Cutoff) ? wrappedPotential.hyperVirial(atoms) : 0;
    }
    
    /**
     * Returns the value of uInt for the wrapped potential.
     */
    public double integral(double rC) {
        return wrappedPotential.integral(rC);
    }

    public double u(double r2) {
        return wrappedPotential.u(r2);
    }

    public double du(double r2) {
        return wrappedPotential.du(r2);
    }

    public double d2u(double r2) {
        return wrappedPotential.d2u(r2);
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
        return new P0Lrc(space, wrappedPotential, this, types);
    }

    /**
     * Inner class that implements the long-range correction for this truncation scheme.
     */
    private static class P0Lrc extends Potential0Lrc implements Potential2Soft {

        private static final long serialVersionUID = 1L;
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
            throw new RuntimeException("Should not be calling gradient on zero-body potential");
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
            integral = -A*Math.pow(rCutoff, D+1)*((Potential2SoftSpherical)truncatedPotential).du(rCutoff*rCutoff) - (D+1)*integral;
            return pairDensity*integral;
        }

        public double integral(double rC) {
            throw new RuntimeException("nope");
        }

        public double u(double r2) {
            throw new RuntimeException("nope");
        }

        public double du(double r2) {
            throw new RuntimeException("nope");
        }

        public double d2u(double r2) {
            throw new RuntimeException("nope");
        }
    }//end of P0lrc
}
