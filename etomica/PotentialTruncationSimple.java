package etomica;       
/**
 * Simple truncation of the potential at an adjustable cutoff separation.
 * The energy is unaffected for separations less than the truncation distance,
 * and it is set to zero beyond this distance.
 *
 * @author David Kofke
 */
public final class PotentialTruncationSimple extends PotentialTruncation {
        
    double rCutoff, r2Cutoff, rCD;
    private /*final*/ double A; //inner class doesn't permit finals
    private /*final*/ int D;
    
    public PotentialTruncationSimple(Potential2 potential) {
        this(potential, Default.POTENTIAL_CUTOFF_FACTOR * Default.ATOM_SIZE);}
        
    public PotentialTruncationSimple(Potential2 potential, double rCutoff) {
        super(potential);
        A = potential.parentSimulation().space().sphereArea(1.0);  //multiplier for differential surface element
        D = potential.parentSimulation().space().D();              //spatial dimension
        setTruncationRadius(rCutoff);
    }
    public boolean isZero(double r2) {return r2 > r2Cutoff;}
        
    public double uTransform(double r2, double untruncatedValue) {
        return (r2 > r2Cutoff) ? 0.0 : untruncatedValue;
    }
    public double duTransform(double r2, double untruncatedValue) {
        return (r2 > r2Cutoff) ? 0.0 : untruncatedValue;
    }
    public double d2uTransform(double r2, double untruncatedValue) {
        return (r2 > r2Cutoff) ? 0.0 : untruncatedValue;
    }
    
    /**
     * Mutator method for the radial cutoff distance.
     */
    public final void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
        rCD = 1.0;
        for(int i=D; i>0; i--) {rCD *= rCD;}  //rC^D
    }
    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {return rCutoff;}
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public etomica.units.Dimension getTruncationRadiusDimension() {return etomica.units.Dimension.LENGTH;}
    
    /**
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(PotentialGroup parent) {
        return new P0Lrc(parent);
    }
    
    /**
     * Inner class that implements the long-range correction for this truncation scheme.
     */
    private class P0Lrc extends Potential0Lrc {
        
        private Phase phase;
        
        public P0Lrc(PotentialGroup parent) {super(parent);}
        public Potential set(Atom a) {return this;}
        public Potential set(Atom a1, Atom a2) {return this;}
        public Potential set(SpeciesMaster s) {return this;}
        
        public double energy() {
            return uCorrection(potential.iterator().size()/phase.volume());
        }
        
        public Potential set(Phase p) {
            phase = p;
            potential.set(p.speciesMaster);
            return this;
        }
        
        /**
         * Long-range correction to the energy.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
        public double uCorrection(double pairDensity) {
            double integral = ((Potential2SoftSpherical)potential).integral(rCutoff);
            return pairDensity*integral;
        }
        
        /**
         * Uses result from integration-by-parts to evaluate integral of
         * r du/dr using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         */
            //not checked carefully
        public double duCorrection(double pairDensity) {
            Potential2SoftSpherical potentialSpherical = (Potential2SoftSpherical)potential;
            double integral = potentialSpherical.integral(rCutoff);
            integral = A*rCD*potentialSpherical.u(rCutoff) - D*integral;//need potential to be spherical to apply here
            return pairDensity*integral;
        }

        /**
         * Uses result from integration-by-parts to evaluate integral of
         * r^2 d2u/dr2 using integral of u.
         * @param pairDensity average pairs-per-volume affected by the potential.
         *
         * Not implemented: always returns zero.
         */
        public double d2uCorrection(double pairDensity) {
            return 0.0;
        }
    }//end of P0lrc
        
}//end of PotentialTruncationSimple
