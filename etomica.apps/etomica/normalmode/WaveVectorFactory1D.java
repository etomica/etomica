package etomica.normalmode;

import java.io.Serializable;

import etomica.api.IVector;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

/**
 * Wave vector factory that returns wave vectors for a 1D system.  
 * These wave vectors are given by 2 Pi m / (N a) = 2 Pi m / L = 2 Pi m rho / N,
 * for m = 1, 2, 3,...,N/2.  If N is odd there will be (N-1)/2 vectors, each with 
 * a coefficient of 1.0 (because none are on the Brillouin-zone boundary) while if 
 * N is even there will be N/2, with the last one having a coefficient of 0.5. Wave
 * vectors corresponding to m <= 0 are not included, although they are technically are among the full
 * set of wave vectors for the system.
 *
 * @author Andrew Schultz and David Kofke
 */
public class WaveVectorFactory1D implements WaveVectorFactory, Serializable {

    public WaveVectorFactory1D(int dim) {
        if(dim != 1) {
            throw new RuntimeException("Must give a box for a 1D system"); 
        }
    }
    
    public void makeWaveVectors(Box box) {

        int nA = box.moleculeCount();
        double L = box.getBoundary().getDimensions().x(0);
        
        int mMax = nA/2;

        waveVectors = new Vector1D[mMax+1];
        coefficients = new double[mMax+1];
        
        for(int m = 0; m<=mMax; m++) {
            waveVectors[m] = new Vector1D(2. * Math.PI * m / L);
            coefficients[m] = 1.0;
        }
        coefficients[0] = 0.5;
        
        if(nA % 2 == 0) {
            coefficients[mMax] = 0.5;
        }

    }
    
    public IVector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getCoefficients() {
        return coefficients;
    }
    
    public static void main(String[] args) {
        int nCells = 6;
        Simulation sim = new Simulation(Space1D.getInstance());
        Box box = new Box(sim, sim.getSpace());
        sim.addBox(box);
        box.setDimensions(new Vector1D(nCells));
        ISpecies species = new SpeciesSpheresMono(sim);
        box.setNMolecules(species, nCells*nCells*nCells);
        
        WaveVectorFactory1D foo = new WaveVectorFactory1D(sim.getSpace().D());
        foo.makeWaveVectors(box);
        IVector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;

    protected IVector[] waveVectors;
    protected double[] coefficients;
}
