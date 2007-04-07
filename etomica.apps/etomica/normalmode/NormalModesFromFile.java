package etomica.normalmode;

import etomica.phase.Phase;

/**
 * Provides normal-mode information as obtained from a file read at
 * construction.
 * 
 * @author kofke
 * 
 */
public class NormalModesFromFile implements NormalModes {

    /**
     * @param filename
     *            Root of file name where normal mode information may be found.
     *            Wavevectors will be taken filename.k, eigenvalues from
     *            filename.val, and eigenvectors from filename.vec
     * @param D spatial dimension of the system (not necessarily coordinate dimension).
     */
    public NormalModesFromFile(String filename, int D) {
        eigenvalues = ArrayReader1D.getFromFile(filename + ".val");
        eigenvectors = ArrayReader2D.getFromFile(filename + ".vec");
        waveVectorFactory = new WaveVectorFactoryFromFile(filename, D);
    }

    public double[][] getEigenvalues(Phase phase) {
        return eigenvalues;
    }

    public double[][][] getEigenvectors(Phase phase) {
        return eigenvectors;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }

    double[][] eigenvalues;
    double[][][] eigenvectors;
    WaveVectorFactory waveVectorFactory;

}
