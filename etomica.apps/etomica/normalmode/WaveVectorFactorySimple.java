package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Wave vector factory that returns wave vectors appropriate for a phase with 
 * a single-atom basis.  The box shape need not be rectangular so long as it
 * matches the primitive's shape.
 *
 * @author Andrew Schultz
 */
public class WaveVectorFactorySimple implements WaveVectorFactory, Serializable {

    public Vector[] makeWaveVectors(Phase phase, Primitive primitive) {
        // If we weren't given wave vectors, determine them from the phase boundary and primitve
        // assume 1-molecule basis and matchup betwen the box and the primitive
    
        double[] d = primitive.getSize();
        int[] numCells = new int[phase.space().D()];
        Vector inverseDim = phase.space().makeVector();
        inverseDim.E(2*Math.PI);
        inverseDim.DE(phase.getBoundary().getDimensions());
        
        for (int i=0; i<phase.space().D(); i++) {
            numCells[i] = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (d[i]));
        }
    
        // 0 to halfSize, round down for halfSize
        // double count here, and then we'll fix it later
        int numWaveVectors = 2*(numCells[0]/2);
        for (int i=0; i<phase.space().D(); i++) {
            // -halfSize to halfSize in the other directions, including 0
            // round down for halfSize
            numWaveVectors *= 2*(numCells[i]/2)+1;
        }
        // exclude <0, 0, 0>
        numWaveVectors--;
        // we exclude negatives of other wave vectors
        numWaveVectors /= 2;
    
        Vector[] waveVectors = new Vector[numWaveVectors];
        int count = 0;
        for (int kx = 0; kx < numCells[0]/2+1; kx++) {
            for (int ky = ((kx==0) ? 1 : -numCells[1]/2); ky < numCells[1]/2+1; ky++) {
                for (int kz = ((kx==0 && ky==0) ? 1 : -numCells[2]/2); kz < numCells[2]/2+1; kz++) {
                    waveVectors[count] = Space.makeVector(new double[]{kx, ky, kz});
                    waveVectors[count].TE(inverseDim);
                    count++;
                }
            }
        }
        return waveVectors;
    }

    private static final long serialVersionUID = 1L;
}
