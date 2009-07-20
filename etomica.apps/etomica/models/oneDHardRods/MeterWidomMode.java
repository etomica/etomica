package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.units.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * 
 * @author cribbin
 *
 */
public class MeterWidomMode extends DataSourceScalar {

    public int nInsert;
    public int affectedWV;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition coordinateDefinition;
    private int coordinateDim;
    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    protected double temperature;
    private double[] uNow, deltaU;
    private double[] waveVectorCoefficients;
    private double wvc;
    
    
    public MeterWidomMode(String string, IPotentialMaster 
            potentialMaster, CoordinateDefinition cd, IBox box){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
    }
    
    
    
    @Override
    public double getDataAsScalar() {
        
        
        return 0;
    }

    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(IVectorMutable[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setAffectedWV(int awv) {
        this.affectedWV = awv;
        wvc = waveVectorCoefficients[awv];
    }
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    public void setWaveVectorCoefficients(double[] waveVectorCoefficients) {
        this.waveVectorCoefficients = waveVectorCoefficients;
    }
    public void setCoordinateDefinition(CoordinateDefinition cd){
        coordinateDefinition = cd;
        coordinateDim = coordinateDefinition.getCoordinateDim();
    }
    
    public void setSpringConstants(double[][] sc){
        omegaSquared = sc;
    }
    
    public void setOmegaSquared(double[][] sc){
        omegaSquared = sc;
    }
}
