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

 * 
 * @author cribbasket
 *
 */
public class MeterNMCBaskets implements IEtomicaDataSource {

    private double eigenVectors[][][];
    private Vector[] waveVectors;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] omegaSquared;
    private int coordinateDim, jump;
    private DataDoubleArray coords;
    
    int nBaskets;
    private int[][] baskets;
    private boolean automaticBaskets;
    
    protected final DataTag tag;
    protected final DataInfoDoubleArray dataInfo;
    private static final long serialVersionUID = 1L;

    public MeterNMCBaskets(CoordinateDefinition coordinateDefinition,
                           Vector[] wv, int nBaskets){
        this(coordinateDefinition, wv, nBaskets, true);
    }
    
    public MeterNMCBaskets(CoordinateDefinition coordinateDefinition,
                           Vector[] wv, int nBaskets, boolean automaticBaskets){
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
        
        if (wv.length % nBaskets != 0 ){
            throw new IllegalStateException("uneven number of wavevectors per basket in MeterMNCBaskets");
        }
        this.nBaskets = nBaskets;
        baskets = new int[nBaskets][wv.length/nBaskets];
        if(automaticBaskets) { automaticBasketsFill(); }
    }
    
    public DataDoubleArray getData(){
        
        for (int basketCount = 0; basketCount < nBaskets; basketCount++){
            for(int incount = 0; incount < baskets[0].length; incount++){
                int wvCount = baskets[basketCount][incount];
            
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
        }
        return coords;
    }
    
    /**
     * sets up the baskets automatically, instead of having the user assign 
     * each wavevector to a basket.
     * 
     */
    private void automaticBasketsFill(){
        int wvPerBasket = baskets[0].length;
        int nWV = baskets.length * baskets[0].length;
        for (int k = 0; k < nWV; k++){
            for (int j = 0; j < nBaskets; j++){
                for (int i = 0; i < wvPerBasket; i++){ 
                    baskets[j][i] = k;
                }
            }
        }
        
    }

    /**
     * Allows the user to choose which wavevectors are in which basket.
     * 
     */
    public void assignToBasket(int whichBasket, int[] theseWVs){
        baskets[whichBasket] = theseWVs;
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


    
}
