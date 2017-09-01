/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.space.Vector;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.normalmode.CoordinateDefinition;
import etomica.units.dimensions.Null;

/**
 * Class which calculates all of the normal mode coordinates for a system.  
 * Data is returned in an array; the first half (entries [0] to 
 * [waveVectors.length*coordinateDim -1]) of which is the real coordinates, 
 * the second half (entries [waveVectors.length*coordinateDim] to 
 * [2*waveVectors.length*coordinateDim]) is the imaginary coordinates.
 * 
 * @author cribbin
 *
 */
public class MeterNormalModeCoordinate implements IEtomicaDataSource {

    public MeterNormalModeCoordinate(CoordinateDefinition coordinateDefinition, Vector[] wv){
        this.coordinateDefinition = coordinateDefinition;
        coordinateDim = this.coordinateDefinition.getCoordinateDim();
        this.waveVectors = wv;
        coords = new DataDoubleArray(waveVectors.length*2*coordinateDim);
        
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        jump = coordinateDim * waveVectors.length;
        
        dataInfo = new DataInfoDoubleArray("Real and Imaginary Normal Mode " +
                "Coordinates", Null.DIMENSION, new int[]{waveVectors.length*2*coordinateDim});
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
            }
            
            for(int j = 0; j < coordinateDim; j++){
                values[j + coordinateDim * wvCount] = realCoord[j];
                values[j + coordinateDim * wvCount + jump] = imagCoord[j];
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

    public int getJump(){
        return jump;
    }


    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }

    private double eigenVectors[][][];
    private Vector[] waveVectors;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] omegaSquared;
    private int coordinateDim, jump;
    private DataDoubleArray coords;

    protected final DataTag tag;
    protected final DataInfoDoubleArray dataInfo;
    private static final long serialVersionUID = 1L;
    
}
