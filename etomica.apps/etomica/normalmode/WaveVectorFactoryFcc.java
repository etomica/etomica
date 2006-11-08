package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
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

    public Vector[] makeWaveVectors(Phase phase, Primitive primitive) {
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
        
        // FCC has 4-atom basis
        int numWaveVectors = 4;
        for (int i=0; i<phase.space().D(); i++) {
            // -halfSize to halfSize in the other directions, including 0
            numWaveVectors *= numCells;
        }
        
        //XXX the given constraints are for FCC
        Vector[] waveVectorsTemp = new Vector[numWaveVectors];
        int count = 0;
        for (int kx = 0; kx <= numCells; kx++) {
            for (int ky = ((kx==0) ? 0 : -numCells); ky <= numCells; ky++) {
                for (int kz = ((kx==0 && ky==0) ? 1 : -numCells); kz <= numCells; kz++) {
                    if (2 * (kx + ky + kz) <= 3 * numCells
                            && 2 * (kx + ky + kz) > -3 * numCells
                            && 2 * (kx + ky - kz) <= 3 * numCells
                            && 2 * (kx + ky - kz) > -3 * numCells
                            && 2 * (kx - ky + kz) <= 3 * numCells
                            && 2 * (kx - ky + kz) > -3 * numCells
                            && 2 * (kx - ky - kz) <= 3 * numCells
                            && 2 * (kx - ky - kz) > -3 * numCells) {
                        waveVectorsTemp[count] = new Vector3D(kx, ky, kz);
                        waveVectorsTemp[count].TE(Math.sqrt(2) * Math.PI / d / numCells);
                        count++;
                    }
                }
            }
        }
        numWaveVectors = count;
        Vector[] waveVectors = new Vector[numWaveVectors];
        System.arraycopy(waveVectorsTemp,0,waveVectors,0,numWaveVectors);

        return waveVectors;
    }
    
    public static void main(String[] args) {
        int nCells = 2;
        Simulation sim = new Simulation(Space3D.getInstance());
        Phase phase = new Phase(sim);
        phase.setDimensions(new Vector3D(nCells, nCells, nCells));
        Species species = new SpeciesSpheresMono(sim);
        phase.getAgent(species).setNMolecules(4*nCells*nCells*nCells);
        Primitive primitive = new PrimitiveCubic(sim.space, 1/Math.sqrt(2));
        
        WaveVectorFactoryFcc foo = new WaveVectorFactoryFcc();
        Vector[] waveVectors = foo.makeWaveVectors(phase, primitive);
        System.out.println("number of wave vectors "+waveVectors.length);
        for (int i=0; i<waveVectors.length; i++) {
            System.out.println(waveVectors[i]);
        }
    }
    
    private static final long serialVersionUID = 1L;
}
