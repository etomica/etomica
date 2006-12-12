package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Wave vector factory that returns wave vectors for a 1D system.
 *
 * @author Andrew Schultz
 */
public class WaveVectorFactory1D implements WaveVectorFactory, Serializable {

    public WaveVectorFactory1D(Primitive primitive) {
        this.primitive = primitive;
        if (primitive.getSpace().D() != 1) {
            throw new IllegalArgumentException("must be 1D");
        }
    }
    
    public void makeWaveVectors(Phase phase) {
        // If we weren't given wave vectors, determine them from the phase boundary and primitve
        // assume 1-molecule basis and matchup betwen the box and the primitive
    
        double[] d = primitive.getSize();
        int[] numCells = new int[phase.space().D()];
        Vector[] reciprocals =  primitive.reciprocal().vectors();;
        Vector[] waveVectorBasis = new Vector[reciprocals.length];
        
        for (int i=0; i<phase.space().D(); i++) {
            waveVectorBasis[i] = phase.space().makeVector();
            waveVectorBasis[i].E(reciprocals[i]);
            numCells[i] = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (d[i]));
            waveVectorBasis[i].TE(1.0/numCells[i]);
        }
    
        int[] kMin = new int[phase.space().D()];
        int[] kMax= new int[phase.space().D()];
        for (int i=0; i<kMax.length; i++) {
            kMin[i] = -(numCells[i]-1)/2;
            kMax[i] = numCells[i]/2;
        }
        int[] waveVectorIndices = new int[2*kMax[0]+1];
        int count = 0;
        // this will find numCells vectors.  Some of them have negatives 
        // within the set others do not.  If its negative is within the set, 
        // exclude the negative, but remember it was there -- they will have 
        // coefficients of '1' while the ones without a negative in the set 
        // will have coefficients of '0.5'.  The ones without a negative have
        // instead a vector which handles the same degree of freedom.
        for (int kx = kMin[0]; kx < kMax[0]+1; kx++) {
            if (kx == 0) continue;
            boolean flip = kx < 0;
            if (flip) {
                if (waveVectorIndices[-kx+kMax[0]] == 0) {
                    // this one was unique
                    count++;
                }
                waveVectorIndices[-kx+kMax[0]]++;
            }
            else {
                if (waveVectorIndices[kx+kMax[0]] == 0) {
                    // this one was unique
                    count++;
                }
                waveVectorIndices[kx+kMax[0]]++;
            }
        }
        waveVectors = new Vector1D[count];
        coefficients = new double[count];
        count = 0;
        for (int kx = -kMax[0]; kx < kMax[0]+1; kx++) {
            if (waveVectorIndices[kx+kMax[0]] > 0) {
                waveVectors[count] = phase.space().makeVector();
                waveVectors[count].Ea1Tv1(kx, waveVectorBasis[0]);
                coefficients[count] = waveVectorIndices[kx+kMax[0]]/2.0;
                count++;
            }
        }
    }
    
    public Vector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getCoefficients() {
        return coefficients;
    }
    
    public static void main(String[] args) {
        int nCells = 6;
        Simulation sim = new Simulation(Space1D.getInstance());
        Phase phase = new Phase(sim);
        phase.setDimensions(new Vector1D(nCells));
        Species species = new SpeciesSpheresMono(sim);
        phase.getAgent(species).setNMolecules(nCells*nCells*nCells);
        Primitive primitive = new PrimitiveCubic(sim.getSpace(), 1);
        
        WaveVectorFactory1D foo = new WaveVectorFactory1D(primitive);
        foo.makeWaveVectors(phase);
        Vector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;
    protected final Primitive primitive;
    protected Vector[] waveVectors;
    protected double[] coefficients;
}
