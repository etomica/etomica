/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.IDataInfo;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Length;

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
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public IDataInfo getDataInfo() {
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
    	return realCoord;
//    	double QCoord = Math.sqrt(realCoord*realCoord + imaginaryCoord*imaginaryCoord);
//    	return realCoord > 0 ? QCoord : -QCoord; 
       
    }

    
    public void setEigenvectors(double[] eigenVectors) {
        this.eigenvectors = eigenVectors;
    }
    
    public void setWaveVector(Vector waveVector){
    	this.waveVector = waveVector;
    } 

    
    
    private static final long serialVersionUID = 1L;
    protected final CoordinateDefinition coordinateDefinition;
    protected double[] eigenvectors;
    protected double[] realT, imaginaryT;
    protected Vector waveVector;
}
