package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.units.Energy;

/**
 * Meter that calculates the Boltzmann-factored harmonic energy of each normal mode for a 
 * configuration given eigenvectors and omegas corresponding to wave vectors.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicSingleEnergy implements DataSource {

    public MeterHarmonicSingleEnergy(CoordinateDefinition coordinateDefinition, NormalModes normalModes) {
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;
        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{0});
        tag = new DataTag();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    

    public Data getData() {
        double[] x = data.getData();
        
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            int coordinateDim = coordinateDefinition.getCoordinateDim();
            for (int i=0; i<coordinateDim; i++) {
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += realT[j] * eigenvectors[iVector][j][i];
                    imaginaryCoord += imaginaryT[j] * eigenvectors[iVector][j][i];
                }
                double normalCoord = (realCoord*realCoord + imaginaryCoord*imaginaryCoord);
                x[iVector*coordinateDim+i] = Math.exp(-0.5 * waveVectorCoefficients[iVector] * 
                        normalCoord * omegaSquared[iVector][i] / temperature);
            }
        }
        return data;
    }
    
    public Phase getPhase() {
        return coordinateDefinition.getPhase();
    }

    public void setPhase(Phase newPhase) {
        coordinateDefinition.setPhase(newPhase);
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        
        normalModes.getWaveVectorFactory().makeWaveVectors(newPhase);
        setWaveVectors(normalModes.getWaveVectorFactory().getWaveVectors(),normalModes.getWaveVectorFactory().getCoefficients());
        setEigenvectors(normalModes.getEigenvectors(newPhase));
        setEigenvalues(normalModes.getEigenvalues(newPhase));

        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{waveVectors.length,coordinateDim});
        data = new DataDoubleArray(new int[]{waveVectors.length,coordinateDim});


        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];
    }
    
    public void setWaveVectors(IVector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    public void setEigenvectors(double[][][] eigenVectors) {
        this.eigenvectors = (double[][][])eigenVectors.clone();
    }
    
    public void setEigenvalues(double[][] eigenvalues) {
        omegaSquared = (double[][])eigenvalues.clone();
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = 1.0/omegaSquared[i][j];
            }
        } 
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    private static final long serialVersionUID = 1L;
    protected final CoordinateDefinition coordinateDefinition;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    private final DataTag tag;
    protected double temperature;
    protected double[] realT, imaginaryT;
    protected IVector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenvectors;
    protected double[][] omegaSquared;
    protected String name;
    protected NormalModes normalModes;
}
