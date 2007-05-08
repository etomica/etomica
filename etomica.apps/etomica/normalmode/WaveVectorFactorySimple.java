package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Wave vector factory that returns wave vectors appropriate for a phase with 
 * a single-atom basis.  The box shape need not be rectangular so long as it
 * matches the primitive's shape.
 *
 * @author Andrew Schultz
 */
public class WaveVectorFactorySimple implements WaveVectorFactory, Serializable {

    public WaveVectorFactorySimple(Primitive primitive) {
        this.primitive = primitive;
    }
    
    public void makeWaveVectors(Phase phase) {
        // If we weren't given wave vectors, determine them from the phase boundary and primitve
        // assume 1-molecule basis and matchup betwen the box and the primitive
    
        double[] d = primitive.getSize();
        int[] numCells = new int[phase.getSpace().D()];
        IVector[] reciprocals =  primitive.reciprocal().vectors();
        IVector[] waveVectorBasis = new IVector[reciprocals.length];
        
        for (int i=0; i<phase.getSpace().D(); i++) {
            waveVectorBasis[i] = phase.getSpace().makeVector();
            waveVectorBasis[i].E(reciprocals[i]);
            numCells[i] = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (d[i]));
            waveVectorBasis[i].TE(1.0/numCells[i]);
        }
    
        int[] kMin = new int[phase.getSpace().D()];
        int[] kMax= new int[phase.getSpace().D()];
        for (int i=0; i<kMax.length; i++) {
            kMin[i] = -(numCells[i]-1)/2;
            kMax[i] = numCells[i]/2;
        }
        int[][][] waveVectorIndices = new int[2*kMax[0]+1][2*kMax[1]+1][2*kMax[2]+1];
        int count = 0;
        int[] idx = new int[3];
        // if we have an odd number of cells, we flip the first non-zero component
        // positive.  If we have an even number of cells, we do so, but ignore any
        // component equal to the max.  And, for even number of cells, if we flip,
        // we re-flip any component that flips to -max. 
        boolean flip2 = numCells[0] % 2 == 0;
        // this will find N-1 vectors.  Some of them have negatives 
        // within the set others do not.  If its negative is within the set, 
        // exclude the negative, but remember it was there -- they will have 
        // coefficients of '1' while the ones without a negative in the set 
        // will have coefficients of '0.5'.
        for (int kx = kMin[0]; kx < kMax[0]+1; kx++) {
            for (int ky = kMin[1]; ky < kMax[1]+1; ky++) {
                for (int kz = kMin[2]; kz < kMax[2]+1; kz++) {
                    if (kx == 0 && ky == 0 && kz == 0) continue;
                    for (int i=0; i<3; i++) {
                        idx[i] = kMax[i];
                    }
                    boolean flip = kx < 0 || (kx == 0 && ky < 0) || (kx == 0 && ky == 0 && kz < 0);
                    if (flip2) {
                        flip = kx < 0 || ((kx == 0 || kx == kMax[0]) && ky < 0) || ((kx == 0 || kx == kMax[0]) && (ky == 0 || ky == kMax[1]) && kz < 0);
                    }
                    if (flip) {
                        idx[0] -= kx;
                        idx[1] -= ky;
                        idx[2] -= kz;
                    }
                    else {
                        idx[0] += kx;
                        idx[1] += ky;
                        idx[2] += kz;
                    }
                    if (flip2 && flip) {
                        for (int i=0; i<3; i++) {
                            if (idx[i] == 0) {
                                idx[i] = 2*kMax[i];
                            }
                        }
                    }

                    if (waveVectorIndices[idx[0]][idx[1]][idx[2]] == 0) {
                        // this one was unique
                        count++;
                    }
                    waveVectorIndices[idx[0]][idx[1]][idx[2]]++;
                }
            }
        }
        waveVectors = new Vector3D[count];
        coefficients = new double[count];
        count = 0;
        for (int kx = -kMax[0]; kx < kMax[0]+1; kx++) {
            for (int ky = -kMax[1]; ky < kMax[1]+1; ky++) {
                for (int kz = -kMax[2]; kz < kMax[2]+1; kz++) {
                    if (waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]] > 0) {
                        waveVectors[count] = phase.getSpace().makeVector();
                        waveVectors[count].Ea1Tv1(kx, waveVectorBasis[0]);
                        waveVectors[count].PEa1Tv1(ky, waveVectorBasis[1]);
                        waveVectors[count].PEa1Tv1(kz, waveVectorBasis[2]);
                        coefficients[count] = waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]]/2.0;
                        count++;
                    }
                }
            }
        }
    }
    
    public IVector[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getCoefficients() {
        return coefficients;
    }
    
    public static void main(String[] args) {
        int nCells = 2;
        Simulation sim = new Simulation(Space3D.getInstance());
        Phase phase = new Phase(sim);
        sim.addPhase(phase);
        phase.setDimensions(new Vector3D(nCells, nCells, nCells));
        Species species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);
        phase.getAgent(species).setNMolecules(nCells*nCells*nCells);
        Primitive primitive = new PrimitiveCubic(sim.getSpace(), 1);
        
        WaveVectorFactorySimple foo = new WaveVectorFactorySimple(primitive);
        foo.makeWaveVectors(phase);
        IVector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;
    protected final Primitive primitive;
    protected IVector[] waveVectors;
    protected double[] coefficients;
}
