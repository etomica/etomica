package etomica.normalmode;

import java.util.Arrays;

import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Length;

/**
 * 
 */
 
public class MeterHarmonicCoordinate extends DataSourceScalar {

    public MeterHarmonicCoordinate(CoordinateDefinition coordinateDefinition) {
    	super("Normal Mode Coordinate", Length.DIMENSION);
    	
        this.coordinateDefinition = coordinateDefinition;        
        realT = new double[coordinateDefinition.getCoordinateDim()];
        imaginaryT = new double[coordinateDefinition.getCoordinateDim()];
        
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    

    public double getDataAsScalar() {
    	
    	int coordinateDim = coordinateDefinition.getCoordinateDim();
    	
    	double realCoord =0;
    	double imaginaryCoord =0;
        
    	coordinateDefinition.calcT(waveVector, realT, imaginaryT);
    	
    	for (int i=0; i<coordinateDim; i++) {
    		realCoord += realT[i] * eigenvectors[i];
    		imaginaryCoord += imaginaryT[i] * eigenvectors[i];
    	}
    	double QCoord = Math.sqrt(realCoord*realCoord + imaginaryCoord*imaginaryCoord);
    	return realCoord > 0 ? QCoord : -QCoord; 
       
    }

    
    public void setEigenvectors(double[] eigenVectors) {
        this.eigenvectors = eigenVectors;
    }
    
    public void setWaveVector(IVector waveVector){
    	this.waveVector = waveVector;
    } 

    
    
    private static final long serialVersionUID = 1L;
    protected final CoordinateDefinition coordinateDefinition;
    protected double[] eigenvectors;
    protected double[] realT, imaginaryT;
    protected IVector waveVector;
}
