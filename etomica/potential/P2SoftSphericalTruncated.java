package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncated extends Potential2SoftSpherical
               implements PotentialTruncated {
    
    public P2SoftSphericalTruncated(Potential2SoftSpherical potential, double truncationRadius) {
        super(potential.getSpace());
        this.potential = potential;
        setTruncationRadius(truncationRadius);
    }
    
    /**
     * Returns the wrapped potential.
     */
    public Potential2SoftSpherical getWrappedPotential() {
        return potential;
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
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
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
    
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public etomica.units.Dimension getTruncationRadiusDimension() {return etomica.units.Length.DIMENSION;}
    
    /**
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(AtomType[] types) {
        return new P0Lrc(space, potential, this, types);
    }
    
    /**
     * Inner class that implements the long-range correction for this truncation scheme.
     */
    private static class P0Lrc extends Potential0Lrc {
        
        private static final long serialVersionUID = 1L;
        private final double A;
        private final int D;
        private Potential2Soft potential;
        
        public P0Lrc(Space space, Potential2Soft truncatedPotential, 
                Potential2Soft potential, AtomType[] types) {
            super(space, types, (Potential2SoftSpherical)truncatedPotential);
            this.potential = potential;
            A = space.sphereArea(1.0);  //multiplier for differential surface element
            D = space.D();              //spatial dimension
        }
 
        public double energy(AtomSet atoms) {
            return uCorrection(nPairs()/box.volume());
        }
        
        public double virial(AtomSet atoms) {
            return duCorrection(nPairs()/box.volume());
        }
        
        public IVector[] gradient(AtomSet atoms) {
            throw new RuntimeException("Should not be calling gradient on zero-body potential");
        }
        
        public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
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
         * r^2 d2u/dr2 using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         *
         * Not implemented: throws RuntimeException.
         */
        public double d2uCorrection(double pairDensity) {
            throw new etomica.exception.MethodNotImplementedException();
        }
    }//end of P0lrc
    
    private static final long serialVersionUID = 1L;
    protected double rCutoff, r2Cutoff;
    protected final Potential2SoftSpherical potential;
}
