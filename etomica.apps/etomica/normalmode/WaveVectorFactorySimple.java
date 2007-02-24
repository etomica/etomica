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
        int[] numCells = new int[phase.space().D()];
        IVector[] reciprocals =  primitive.reciprocal().vectors();;
        IVector[] waveVectorBasis = new IVector[reciprocals.length];
        
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
        int[][][] waveVectorIndices = new int[2*kMax[0]+1][2*kMax[1]+1][2*kMax[2]+1];
        int count = 0;
        // this will find 4(numCells)^3 vectors.  Some of them have negatives 
        // within the set others do not.  If its negative is within the set, 
        // exclude the negative, but remember it was there -- they will have 
        // coefficients of '1' while the ones without a negative in the set 
        // will have coefficients of '0.5'.  The ones without a negative have
        // instead a vector which handles the same degree of freedom.
        for (int kx = kMin[0]; kx < kMax[0]+1; kx++) {
            for (int ky = kMin[1]; ky < kMax[1]+1; ky++) {
                for (int kz = kMin[2]; kz < kMax[2]+1; kz++) {
                    if (kx == 0 && ky == 0 && kz == 0) continue;
                    boolean flip = kx < 0 || (kx == 0 && ky < 0) || (kx == 0 && ky == 0 && kz < 0);
                    if (flip) {
                        if (waveVectorIndices[-kx+kMax[0]][-ky+kMax[1]][-kz+kMax[2]] == 0) {
                            // this one was unique
                            count++;
                        }
                        waveVectorIndices[-kx+kMax[0]][-ky+kMax[1]][-kz+kMax[2]]++;
                    }
                    else {
                        if (waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]] == 0) {
                            // this one was unique
                            count++;
                        }
                        waveVectorIndices[kx+kMax[0]][ky+kMax[1]][kz+kMax[2]]++;
                    }
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
                        waveVectors[count] = phase.space().makeVector();
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
        phase.setDimensions(new Vector3D(nCells, nCells, nCells));
        Species species = new SpeciesSpheresMono(sim);
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
