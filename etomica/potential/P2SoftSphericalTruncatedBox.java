package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncatedBox extends Potential2SoftSpherical
               implements PotentialTruncated {
    
    public P2SoftSphericalTruncatedBox(Potential2SoftSpherical potential) {
        super(potential.getSpace());
        this.potential = potential;
        setTruncationFactor(0.49);
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

    public void setBox(IBox newBox) {
        super.setBox(newBox);
        IVector dim = newBox.getBoundary().getDimensions();
        double minL = dim.x(0);
        for (int i=1; i<dim.getD();  i++) {
            if (minL > dim.x(i)) {
                minL = dim.x(i);
            }
        }
        setTruncationRadius(minL*factor);
    }
    
    public void setTruncationFactor(double newFactor) {
        if (newFactor <= 0 || newFactor > 0.5) {
            throw new RuntimeException("Factor must be positive and must not be greater than 0.5");
        }
        factor = newFactor;
    }
    
    public double getTruncationFactor() {
        return factor;
    }
    
    /**
     * Mutator method for the radial cutoff distance.
     */
    protected void setTruncationRadius(double rCut) {
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
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(IAtomType[] types) {
        return new P0Lrc(space, potential, this, types);
    }
    
    /**
     * Inner class that implements the long-range correction for this truncation scheme.
     */
    private static class P0Lrc extends Potential0Lrc {
        
        private static final long serialVersionUID = 1L;
        private final double A;
        private final int D;
        private P2SoftSphericalTruncatedBox potential;
        
        public P0Lrc(Space space, Potential2SoftSpherical truncatedPotential, 
                P2SoftSphericalTruncatedBox potential, IAtomType[] types) {
            super(space, types, truncatedPotential);
            this.potential = potential;
            A = space.sphereArea(1.0);  //multiplier for differential surface element
            D = space.D();              //spatial dimension
        }
 
        public double energy(IAtomSet atoms) {
            return uCorrection(nPairs()/box.volume());
        }
        
        public double virial(IAtomSet atoms) {
            return duCorrection(nPairs()/box.volume());
        }
        
        public IVector[] gradient(IAtomSet atoms) {
            throw new RuntimeException("Should not be calling gradient on zero-body potential");
        }
        
        public IVector[] gradient(IAtomSet atoms, Tensor pressureTensor) {
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
            double integral = ((Potential2SoftSpherical)truncatedPotential).integral(rCutoff);
            return pairDensity*integral;
        }
        
        /**
         * Uses result from integration-by-parts to evaluate integral of
         * r du/dr using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
        public double duCorrection(double pairDensity) {
            double rCutoff = potential.getRange();
            double integral = ((Potential2SoftSpherical)truncatedPotential).integral(rCutoff);
            //need potential to be spherical to apply here
            integral = -A*space.powerD(rCutoff)*((Potential2SoftSpherical)truncatedPotential).u(rCutoff*rCutoff) - D*integral;
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
    protected double factor;
}
