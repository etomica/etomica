package etomica;       
/**
 * Simple truncation of the potential as an adjustable cutoff separation.
 * The energy is unaffected for separations less than the truncation distance,
 * and it is set to zero beyond this distance.
 *
 * @author David Kofke
 */
public final class PotentialTruncationSimple extends PotentialTruncation {
        
    double rCutoff, r2Cutoff, rCD;
    private final double A;
    private final int D;
    private Potential2SoftSpherical potentialSpherical;
        
    public PotentialTruncationSimple(Potential2SoftSpherical potential) {
        this(potential, Default.POTENTIAL_CUTOFF_FACTOR * Default.ATOM_SIZE);}
        
    public PotentialTruncationSimple(Potential2SoftSpherical potential, double rCutoff) {
        super(potential);
        A = potential.parentSimulation().space().sphereArea(1.0);  //multiplier for differential surface element
        D = potential.parentSimulation().space().D();                 //spatial dimension
        setCutoff(rCutoff);
        potentialSpherical = potential;
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
    
    public double uCorrection(double pairDensity) {
        double integral = parentPotential.integral(rCutoff);
        return pairDensity*integral;
    }
    /**
     * Uses result from integration-by-parts to evaluate integral of
     * r du/dr using integral of u.
     */
        //not checked carefully
    public double duCorrection(double pairDensity) {
        double integral = parentPotential.integral(rCutoff);
        integral = A*rCD*potentialSpherical.u(rCutoff) - D*integral;//need potential to be spherical to apply here
        return pairDensity*integral;
    }

    /**
     * Uses result from integration-by-parts to evaluate integral of
     * r^2 d2u/dr2 using integral of u.
     *
     * Not implemented: always returns zero.
     */
    public double d2uCorrection(double pairDensity) {
        return 0.0;
    }
        
    public final void setCutoff(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
        rCD = 1.0;
        for(int i=D; i>0; i--) {rCD *= rCD;}  //rC^D
    }
    public double getCutoff() {return rCutoff;}
    public etomica.units.Dimension getCutoffDimension() {return etomica.units.Dimension.LENGTH;}
        
}//end of PotentialTruncationSimple
