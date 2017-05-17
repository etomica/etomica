/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;

/**
 *
 *
 * @author Tai Boon Tan
 */
public class NormalModes3D implements NormalModes {

    /**
     */
    public NormalModes3D(Space _space, Primitive primitive, Basis basis) {
    	
    	this.space = _space;
    	this.primitive = primitive;
    	this.basis = basis;
    	waveVectorFactory = new WaveVectorFactorySimple(primitive, space);
        harmonicFudge = 1;
        
    }

    public double[][] getOmegaSquared() {
    	eigenvalues = ArrayReader1D.getFromFile("DB_FCC_n12_N"+getNCellNum()+".val");
    	double[][] omega2 = new double[eigenvalues.length][eigenvalues[0].length];
    	for (int i=0; i<omega2.length; i++) {
    		for (int j=0; j<omega2[i].length; j++) {
    			omega2[i][j] = 1.0/eigenvalues[i][j];
               }
    	}

        return omega2;
    }

    public double[][][] getEigenvectors() {
    	eigenvectors = ArrayReader2D.getFromFile("DB_FCC_n12_N"+getNCellNum()+".vec");
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
    
    public int getNCellNum() {
		return nCellNum;
	}

	public void setNCellNum(int cellNum) {
		nCellNum = cellNum;
	}

	public static void main(String[] args){

   		
    }
    
    protected WaveVectorFactory waveVectorFactory;
    protected double harmonicFudge;
    protected double temperature;
    protected Primitive primitive;
    protected Space space;
    protected Basis basis;
    protected double[][][] eigenvectors;
    protected double[][] eigenvalues;
    protected int nCellNum;
}
