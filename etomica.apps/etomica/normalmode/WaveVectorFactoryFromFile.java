package etomica.normalmode;

import etomica.api.IVector;
import etomica.box.Box;
import etomica.space.Space;

/**
 * Obtains wave vectors and coefficients from a file.
 * 
 * @author kofke
 * 
 */
public class WaveVectorFactoryFromFile implements WaveVectorFactory {

    /**
     * @param filename
     *            such that filename.k holds the wave vectors
     * @param D
     *            spatial dimension (not necessarily coordinate dimension)
     */
    public WaveVectorFactoryFromFile(String filename, int D) {
        // read and process wave vectors
        double[][] waveVectorsAndCoefficients = ArrayReader1D
                .getFromFile(filename + ".k");
        waveVectors = new IVector[waveVectorsAndCoefficients.length];
        coefficients = new double[waveVectors.length];
        double[] justWaveVector = new double[D];
        for (int i = 0; i < waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            for (int j = 0; j < D; j++) {
                justWaveVector[j] = waveVectorsAndCoefficients[i][j + 1];
            }
            waveVectors[i] = Space.makeVector(justWaveVector);
        }

    }

    public double[] getCoefficients() {
        return coefficients;
    }

    public IVector[] getWaveVectors() {
        return waveVectors;
    }

    public void makeWaveVectors(Box box) {
        // nothing to do
    }

    private IVector[] waveVectors;
    private double[] coefficients;
}
