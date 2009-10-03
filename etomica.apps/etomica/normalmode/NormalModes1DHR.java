package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesSpheresMono;

/**
 * Normal-mode quantities for a 1-dimensional system of hard rods.  Frequencies are defined
 * in terms of pair correlations (S matrix), which is given using an analytic solution.
 */
public class NormalModes1DHR implements NormalModes {
    
    public NormalModes1DHR(int D) {
        harmonicFudge = 1;
        this.dim = D;
        waveVectorFactory = new WaveVectorFactory1D();
    }
    
    public double[][] getOmegaSquared(IBox box) {
        if(dim != 1) {
            throw new RuntimeException("Must give a box for a 1D system"); 
        }
        int nA = box.getMoleculeList().getMoleculeCount();
        double L = box.getBoundary().getBoxSize().getX(0);
        int mMax = nA/2;
        double[][] omega2= new double[mMax+1][1];
        omega2[0][0] = Double.POSITIVE_INFINITY;
        for(int m=1; m<=mMax; m++) {
            omega2[m][0] = temperature/(harmonicFudge*S1DHR(m, L, nA));
        }
        return omega2;
    }

    public double[][][] getEigenvectors(IBox box) {
        if(dim != 1) {
            throw new RuntimeException("Must give a box for a 1D system"); 
        }
        int nA = box.getMoleculeList().getMoleculeCount();
        int mMax = nA/2;
        double[][][] eVecs = new double[mMax+1][1][1];
        for(int m=0; m<=mMax; m++) {
            eVecs[m][0][0] = 1.0;  
        }
        return eVecs;
    }

    public void setHarmonicFudge(double newHarmonicFudge) {
        harmonicFudge = newHarmonicFudge;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    protected double harmonicFudge;
    protected double temperature;
    private final WaveVectorFactory1D waveVectorFactory;
    private final int dim;

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

    public static void main(String[] args) {
        int N = 8;
        NormalModes1DHR nm = new NormalModes1DHR(1);
        nm.setTemperature(1);
        ISpace space = Space1D.getInstance();
        Box box = new Box(space);
        Simulation sim = new Simulation(space);
        sim.addBox(box);
        ISpecies species = new SpeciesSpheresMono(sim,space);
        sim.getSpeciesManager().addSpecies(species);
        box.setNMolecules(species, N);
        box.getBoundary().setBoxSize(new Vector1D(2*N));
        double[][] omega2 = nm.getOmegaSquared(box);
        WaveVectorFactory1D wvf = new WaveVectorFactory1D();
        wvf.makeWaveVectors(box);
        for (int i=0; i<omega2.length; i++) {
            System.out.println(wvf.getWaveVectors()[i]+" "+omega2[i][0]);
        }
    }
}
