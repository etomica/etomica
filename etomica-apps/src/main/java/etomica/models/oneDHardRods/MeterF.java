/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.normalmode.CoordinateDefinition;
import etomica.units.dimensions.Null;

public class MeterF extends DataSourceScalar {

    CoordinateDefinition coordinateDefinition;
    Vector[] waveVectors;
    double[] realT, imagT;
    int coordinateDim, convertedWV;
    double[][][] eigenVectors;
    
    public MeterF(CoordinateDefinition cd){
        super("MeterF", Null.DIMENSION);
        coordinateDefinition = cd;
        coordinateDim = coordinateDefinition.getCoordinateDim();
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        
    }
    
    public double getDataAsScalar(){
        coordinateDefinition.calcT(waveVectors[convertedWV], realT, imagT);
        double realCoord = 0.0, imagCoord = 0.0;
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            for(int j = 0; j < coordinateDim; j++){
                realCoord += eigenVectors[convertedWV][i][j] * realT[j];
                imagCoord += eigenVectors[convertedWV][i][j] * imagT[j];
            }
        }
        
        return realCoord*realCoord + imagCoord*imagCoord;
    }

    public void setWaveVectors(Vector[] waveVectors) {
        this.waveVectors = waveVectors;
    }

    public void setConvertedWV(int convertedWV) {
        this.convertedWV = convertedWV;
    }

    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    
    
    
}
