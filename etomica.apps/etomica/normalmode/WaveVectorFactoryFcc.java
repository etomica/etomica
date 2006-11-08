package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.PrimitiveFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * WaveVectorFactory implementation that returns wave vectors appropriate for 
 * a cubic FCC lattice.
 * 
 * @author Andrew Schultz
 */
public class WaveVectorFactoryFcc implements WaveVectorFactory, Serializable {

    public WaveVectorFactoryFcc(PrimitiveFcc primitive) {
        this.primitive = primitive;
    }
    
    public void makeWaveVectors(Phase phase) {
        int numCells = 0;
        double d = -1;
        for (int i=0; i<phase.space().D(); i++) {
            //XXX divide by sqrt(2) for FCC
            int n = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (primitive.getSize()[i]*Math.sqrt(2)));
            if (i>0 && n != numCells) {
                throw new RuntimeException("Things would be so much happier if you would just use the same number of cells in each direction.");
            }
            numCells = n;
            d = primitive.getSize()[i];
        }
        
        int[][][] waveVectorIndices = new int[2*numCells+1][2*numCells+1][2*numCells+1];
        int count = 0;
        // this will find 4(numCells)^3 vectors.  Some of them have negatives 
        // within the set others do not.  If its negative is within the set, 
        // exclude the negative, but remember it was there -- they will have 
        // coefficients of '1' while the ones without a negative in the set 
        // will have coefficients of '0.5'.  The ones without a negative have
        // instead a vector which handles the same degree of freedom.
        for (int kx = -numCells+1; kx <= numCells; kx++) {
            for (int ky = -numCells+1; ky <= numCells; ky++) {
                for (int kz = -numCells+1; kz <= numCells; kz++) {
                    if (kx == 0 && ky == 0 && kz == 0) continue;
                    if (2 * (Math.abs(kx) + Math.abs(ky) + Math.abs(kz)) <= 3 * numCells
                            && 2 * (kx + ky + kz) > -3 * numCells
                            && 2 * (kx + ky - kz) <= 3 * numCells
                            && 2 * (kx + ky - kz) > -3 * numCells
                            && 2 * (kx - ky + kz) <= 3 * numCells
                            && 2 * (kx - ky + kz) > -3 * numCells
                            && 2 * (kx - ky - kz) <= 3 * numCells
                            && 2 * (kx - ky - kz) > -3 * numCells) {
                        
                        boolean flip = kx < 0 || (kx == 0 && ky < 0) || (kx == 0 && ky == 0 && kz < 0);
                        if (flip) {
                            if (waveVectorIndices[-kx+numCells][-ky+numCells][-kz+numCells] == 0) {
                                // this one was unique
                                count++;
                            }
                            waveVectorIndices[-kx+numCells][-ky+numCells][-kz+numCells]++;
                        }
                        else {
                            if (waveVectorIndices[kx+numCells][ky+numCells][kz+numCells] == 0) {
                                // this one was unique
                                count++;
                            }
                            waveVectorIndices[kx+numCells][ky+numCells][kz+numCells]++;
                        }
                    }
                }
            }
        }
        waveVectors = new Vector3D[count];
        coefficients = new double[count];
        count = 0;
        for (int kx = -numCells; kx <= numCells; kx++) {
            for (int ky = -numCells; ky <= numCells; ky++) {
                for (int kz = -numCells; kz <= numCells; kz++) {
                    if (waveVectorIndices[kx+numCells][ky+numCells][kz+numCells] > 0) {
                        waveVectors[count] = new Vector3D(kx, ky, kz);
                        waveVectors[count].TE(Math.sqrt(2) * Math.PI / d / numCells);
                        coefficients[count] = waveVectorIndices[kx+numCells][ky+numCells][kz+numCells]/2.0;
                        count++;
                    }
                }
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
        int nCells = 3;
        Simulation sim = new Simulation(Space3D.getInstance());
        Phase phase = new Phase(sim);
        phase.setDimensions(new Vector3D(nCells, nCells, nCells));
        Species species = new SpeciesSpheresMono(sim);
        phase.getAgent(species).setNMolecules(4*nCells*nCells*nCells);
        PrimitiveFcc primitive = new PrimitiveFcc(sim.space, 1/Math.sqrt(2));
        
        WaveVectorFactoryFcc foo = new WaveVectorFactoryFcc(primitive);
        foo.makeWaveVectors(phase);
        Vector[] waveVectors = foo.getWaveVectors();
        double[] coefficients = foo.getCoefficients();
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(coefficients[i]+" "+waveVectors[i]);
        }
    }

    private static final long serialVersionUID = 1L;
    protected final PrimitiveFcc primitive;
    protected Vector3D[] waveVectors;
    protected double[] coefficients;
}
