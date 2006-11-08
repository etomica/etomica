package etomica.normalmode;

import java.io.Serializable;

import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
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
        Vector[] reciprocals =  primitive.reciprocal().vectors();;
        Vector[] waveVectorBasis = new Vector[reciprocals.length];
        
        for (int i=0; i<phase.space().D(); i++) {
            waveVectorBasis[i] = phase.space().makeVector();
            waveVectorBasis[i].E(reciprocals[i]);
            numCells[i] = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (d[i]));
            waveVectorBasis[i].TE(1.0/numCells[i]);
        }
    
        // 0 to halfSize, round down for halfSize
        // double count here, and then we'll fix it later
        int numWaveVectors = 2*(numCells[0]/2)+1;
        for (int i=1; i<phase.space().D(); i++) {
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
        Vector waveVectorX = phase.space().makeVector();
        Vector waveVectorXY = phase.space().makeVector();
        for (int kx = 0; kx < numCells[0]/2+1; kx++) {
            waveVectorX.Ea1Tv1(kx, waveVectorBasis[0]);
            for (int ky = ((kx==0) ? 0 : -numCells[1]/2); ky < numCells[1]/2+1; ky++) {
                waveVectorXY.Ev1Pa1Tv2(waveVectorX, ky, waveVectorBasis[1]);
                for (int kz = ((kx==0 && ky==0) ? 1 : -numCells[2]/2); kz < numCells[2]/2+1; kz++) {
                    waveVectors[count] = phase.space().makeVector();
                    waveVectors[count].Ev1Pa1Tv2(waveVectorXY, kz, waveVectorBasis[2]);
                    count++;
                }
            }
        }
        return waveVectors;
    }

    private static final long serialVersionUID = 1L;
}
