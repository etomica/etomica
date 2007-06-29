package etomica.normalmode;

import etomica.box.Box;

/**
 * Provides normal-mode information as obtained from a file read at
 * construction.  File is expected to contain information about eigenvectors
 * and eigenvalues of pair correlations.  Reciprocal of eigenvalues gives frequencies
 * defined by this implementation.
 * 
 * @author kofke
 * 
 */
public class NormalModesFromFile implements NormalModes {

    /**
     * @param filename
     *            Root of file name where normal mode information may be found.
     *            Wavevectors will be taken filename.k, eigenvalues from
     *            filename.val, and eigenvectors from filename.vec.
     * @param D spatial dimension of the system (not necessarily coordinate dimension).
     */
    public NormalModesFromFile(String filename, int D) {
        double[][] eigenvalues = ArrayReader1D.getFromFile(filename + ".val");
        omega2 = new double[eigenvalues.length][eigenvalues[0].length];
        for (int i=0; i<omega2.length; i++) {
            for (int j=0; j<omega2[i].length; j++) {
                omega2[i][j] = 1.0/eigenvalues[i][j];
            }
        }
        eigenvectors = ArrayReader2D.getFromFile(filename + ".vec");
        waveVectorFactory = new WaveVectorFactoryFromFile(filename, D);
        harmonicFudge = 1;
        temperature = 1;
    }

    public double[][] getOmegaSquared(Box box) {
        double[][] fudgedOmega2 = new double[omega2.length][omega2[0].length];
        for (int i=0; i<fudgedOmega2.length; i++) {
            for (int j=0; j<fudgedOmega2[i].length; j++) {
                fudgedOmega2[i][j] = temperature*omega2[i][j]/harmonicFudge;
            }
        }
        return fudgedOmega2;
    }

    public double[][][] getEigenvectors(Box box) {
        return eigenvectors;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setHarmonicFudge(double newHarmonicFudge) {
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    double[][] omega2;
    double[][][] eigenvectors;
    WaveVectorFactory waveVectorFactory;
    double harmonicFudge;
    double temperature;
}
