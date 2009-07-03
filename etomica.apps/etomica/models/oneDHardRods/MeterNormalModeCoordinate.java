package etomica.models.oneDHardRods;

import etomica.api.IVectorMutable;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.normalmode.CoordinateDefinition;
import etomica.units.Null;

public class MeterNormalModeCoordinate implements IEtomicaDataSource {

    public MeterNormalModeCoordinate(CoordinateDefinition coordinateDefinition, IVectorMutable[] wv){
        this.coordinateDefinition = coordinateDefinition;
        coordinateDim = this.coordinateDefinition.getCoordinateDim();
        this.waveVectors = wv;
        coords = new DataDoubleArray(waveVectors.length*2*coordinateDim);
        
        
        
        dataInfo = new DataInfoDoubleArray("Real and Imaginary Normal Mode " +
                "Coordinates", Null.DIMENSION, new int[]{2*coordinateDim});
        tag = new DataTag();
    }
    
    public DataDoubleArray getData(){
        for(int wvCount = 0; wvCount < waveVectors.length; wvCount++){
            coordinateDefinition.calcT(waveVectors[wvCount], realT, imagT);
            
            double[] realCoord = new double[coordinateDim];
            double[] imagCoord = new double[coordinateDim];
            double[] values = coords.getData();
            for(int j = 0; j < coordinateDim; j++){
                realCoord[j] = 0.0;
                imagCoord[j] = 0.0;
            }
            
            for(int i = 0; i < coordinateDim; i++){
                if(Double.isInfinite(omegaSquared[wvCount][i])){
                    continue;
                }
                for(int j = 0; j < coordinateDim; j++){
                    realCoord[i] += eigenVectors[wvCount][i][j] * realT[j];
                    imagCoord[i] += eigenVectors[wvCount][i][j] * imagT[j];
                }
                System.out.println(realCoord[i]);
            }
            
            for(int j = 0; j < coordinateDim*2; j++){
                values[2*wvCount + j] = realCoord[j];
                values[2*wvCount+j+coordinateDim] = imagCoord[j];
            }
        }
        
        return coords;
    }
    

    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }

    public void setOmegaSquared(double[][] omegaSquared) {
        this.omegaSquared = omegaSquared;
    }




    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }

    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] omegaSquared;
    private int coordinateDim;
    private DataDoubleArray coords;

    protected final DataTag tag;
    protected final DataInfoDoubleArray dataInfo;
    private static final long serialVersionUID = 1L;
    
}
