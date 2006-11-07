package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

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
            for (int ky = ((kx==0) ? 1 : -numCells + 1); ky <= numCells; ky++) {
                for (int kz = ((kx==0 && ky==0) ? 1 : -numCells + 1); kz <= numCells; kz++) {
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

    private static final long serialVersionUID = 1L;
}
