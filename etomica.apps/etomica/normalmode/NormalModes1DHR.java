package etomica.normalmode;

import etomica.phase.Phase;

/**
 * Normal-mode quantities for a 1-dimensional system of hard rods.  S matrix
 * is given using an analytic solution.
 */
public class NormalModes1DHR implements NormalModes {
    
    public double[][] getEigenvalues(Phase phase) {
        if(phase.getSpace().D() != 1) {
            throw new RuntimeException("Must give a phase for a 1D system"); 
        }
        int nA = phase.moleculeCount();
        double L = phase.getBoundary().getDimensions().x(0);
        int mMax = nA/2;
        double[][] eVals = new double[mMax][1];
        for(int m=1; m<=mMax; m++) {
            eVals[m-1][0] = harmonicFudge*S1DHR(m, L, nA);
        }
        return eVals;
    }

    public double[][][] getEigenvectors(Phase phase) {
        if(phase.getSpace().D() != 1) {
            throw new RuntimeException("Must give a phase for a 1D system"); 
        }
        int nA = phase.moleculeCount();
        int mMax = nA/2;
        double[][][] eVecs = new double[mMax][1][1];
        for(int m=1; m<=mMax; m++) {
            eVecs[m-1][0][0] = 1.0;  
        }
        return eVecs;
    }

    public void setHarmonicFudge(double newHarmonicFudge) {
        harmonicFudge = newHarmonicFudge;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    protected double harmonicFudge;
    private final WaveVectorFactory1D waveVectorFactory = new WaveVectorFactory1D();

    /**
     * Returns the analytical result for the S matrix (just a scalar in this
     * case) for the 1D HR system with PBC
     * 
     * @param m
     *            specification of the wave vector, such that k = 2 Pi m rho/N
     * @param L
     *            lenth of the system
     * @param N
     *            number of rods
     */
    public static double S1DHR(int m, double L, int N) {
        if (m == 0) {
            return (L-N) * (L-N) * (N+1) / 12.0 / N;
        }
        double csc = 1 / Math.sin(m * Math.PI / N);
        return (L-N) * (L-N)
                * (N*(N+2) + Math.cos(2*m*(1+1./N)*Math.PI) - (N+1)*(N+1)*Math.cos(2*m*Math.PI/N))
                * (csc * csc * csc * csc) / (8.*N*N*(2 + 3*N + N*N));
    }

}
