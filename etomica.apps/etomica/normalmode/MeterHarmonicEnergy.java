package etomica.normalmode;

import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.units.Energy;

/**
 * Meter that calculates the harmonic energy of a configuration given
 * eigenvectors and omegas corresponding to wave vectors.
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergy extends DataSourceScalar {

    public MeterHarmonicEnergy(CoordinateDefinition coordinateDefinition) {
        super("Harmonic Energy", Energy.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public double getDataAsScalar() {
        double energySum = 0;
        
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            int coordinateDim = coordinateDefinition.getCoordinateDim();
            for (int i=0; i<coordinateDim; i++) {
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[iVector][j][i] * realT[j];
                    imaginaryCoord += eigenVectors[iVector][j][i] * imaginaryT[j];
                }
                double normalCoord = realCoord*realCoord + imaginaryCoord*imaginaryCoord;
                energySum += waveVectorCoefficients[iVector] * normalCoord * omegaSquared[iVector][i];
            }
        }
        return 0.5*energySum;
    }

    public Phase getPhase() {
        return coordinateDefinition.getPhase();
    }

    public void setPhase(Phase newPhase) {
        coordinateDefinition.setPhase(newPhase);

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];

    }
    
    public void setWaveVectors(IVector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    public void setEigenvectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] newOmegaSquared) {
        omegaSquared = newOmegaSquared;
    }
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected double[] realT, imaginaryT;
    protected IVector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenVectors;
    protected double[][] omegaSquared;
}
