/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.space.Boundary;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;

/**
 *
 * 
 */
public class NormalModes2D implements NormalModes {

    /**
     */
    public NormalModes2D(Space _space, Boundary boundary, Primitive primitive, Basis basis) {
    	
    	this.space = _space;
    	this.boundary = boundary;
    	this.primitive = primitive;
    	this.basis = basis;
    	waveVectorFactory = new WaveVectorFactory2D(primitive, space);
        harmonicFudge = 1;
        
    }

    public double[][] getOmegaSquared() {
    	
    	double[] d = primitive.getSize();
    	int[] numCells = new int[space.D()];
    	
    	for(int i=0; i<space.D(); i++){
    		numCells[i] = (int)Math.round(boundary.getBoxSize().getX(i)/(d[i]));
    		
    	}
    	
    	int numWV = waveVectorFactory.getWaveVectors().length;
    	int numEigen = basis.getScaledCoordinates().length*space.D();
    	
    	double[][] omega2 = new double[numWV][numEigen];
    	
    	/*
    	 *  set all the eigenvalues to 0.001
    	 *  omega2 = 1/eigenvalue
    	 */
    	
        for (int i=0; i<omega2.length; i++) {
            for (int j=0; j<omega2[i].length; j++) {
            	if (i==0 && j<1){
            		omega2[i][j] = Double.POSITIVE_INFINITY;
            		
            	} else {
            		omega2[i][j] = 1000*temperature;
            		
            	}
            }
        }
        return omega2;
    }

    
    
    public double[][][] getEigenvectors() {
    	
    	double[] d = primitive.getSize();
    	int[] numCells = new int[space.D()];
    	
    	for(int i=0; i<space.D(); i++){
    		numCells[i] = (int)Math.round(boundary.getBoxSize().getX(i)/(d[i]));
    		
    	}
    	
    	int numWV = waveVectorFactory.getWaveVectors().length;
    	int numEigen = basis.getScaledCoordinates().length*space.D();
    	
    	eigenvectors = new double[numWV][numEigen][numEigen];
    	
    	for(int i=0; i<numWV; i++){
    		for(int j=0; j<numEigen; j++){
    			for(int k=0; k<numEigen; k++){
    				if(j==k){
    					eigenvectors[i][j][k] = 0.8;
    					
    				} else {
    					eigenvectors[i][j][k] = 0.0;
    					
    				}
    			}
    		}
    	}
    	
        return eigenvectors;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setHarmonicFudge(double newHarmonicFudge) {
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }


    protected double[][][] eigenvectors;
    protected WaveVectorFactory waveVectorFactory;
    protected double harmonicFudge;
    protected double temperature;
    protected Primitive primitive;
    protected Space space;
    protected Basis basis;
    protected Boundary boundary;

}
