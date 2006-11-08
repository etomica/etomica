package etomica.normalmode;

import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.space.Vector;

/**
 * Interface for a class the returns the appropriate wave vectors for a phase
 * and primitive
 *
 * @author Andrew Schultz
 */
public interface WaveVectorFactory {

    /**
     * Returns an array of wave vectors appropraite for the given phase and 
     * primitive.  The wave vectors will not include the 0 vector, any two 
     * vectors which are opposites or any wave vector with a frequency that
     * exceeds half the nearest-neighbor frequency based on the primitve.
     */
    public void makeWaveVectors(Phase phase, Primitive primitive);
    
    public Vector[] getWaveVectors();
    
    public double[] getCoefficients();

}