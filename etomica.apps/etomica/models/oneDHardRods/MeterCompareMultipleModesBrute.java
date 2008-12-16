package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialMaster;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.Null;

public class MeterCompareMultipleModesBrute extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    
    private double eigenVectors[][][];
    private IVector[] waveVectors;
    int[] comparedWVs;
    protected double temperature;
    private double[] waveVectorCoefficients;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double energyHardRod, energyHarmonic;
    
    private static final long serialVersionUID = 1L;
    
    
    
    public MeterCompareMultipleModesBrute(IPotentialMaster potentialMaster, 
            CoordinateDefinition cd, IBox box){
        this("meterCompareMultipleModes", potentialMaster, cd, box);
    }
    
    public MeterCompareMultipleModesBrute(String string, IPotentialMaster
            potentialMaster, CoordinateDefinition cd, IBox box){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
    }
    
    
    
    public double getDataAsScalar(){
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        energyHardRod = 0.0;
        energyHarmonic = 0.0;
        
        
        
        
        
        
        
        
        
        return 0.0;
    }
    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(IVector[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setComparedWVs(int[] cwvs) {
        this.comparedWVs = cwvs;
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
