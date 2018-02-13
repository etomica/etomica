/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.molecule.iterator.MoleculeIteratorAllMolecules;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

/**
 * 
 * Integrator that is used for the harmonic osciallating system
 * 
 * 
 * @author Tai Boon Tan
 */
public class IntegratorHarmonic extends IntegratorMD {

    public IntegratorHarmonic(IRandom random, double timeStep, double temperature, Space _space, Box box) {
        super(null,random, timeStep, temperature, _space, box);
        iterator = new MoleculeIteratorAllMolecules();
        // make IntergratorMD happy.
        meterKE = new DataSourceScalar("", Null.DIMENSION) {
            public double getDataAsScalar() {return 0;}
        };
    }
    
    public void doThermostat() {
        // do nothing -- we don't care about velocities
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setOmegaSquared(double[][] omega2, double[] coeff) {
        omega = new double[omega2.length][omega2[0].length];
   
       	for (int i=0; i<omega.length; i++) {
       		for (int j=0; j<omega[i].length; j++) {
       			omega[i][j] = Math.sqrt(omega2[i][j]);
       		} 
       	}
    }
    
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public Vector[] getWaveVectors(){
    	return waveVectors;
    }
    
    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        waveVectorCoefficients = newWaveVectorCoefficients;
    }
    
    public void setEigenVectors(double[][][] newEigenVectors) {
   		eigenVectors = newEigenVectors;
    }
    
    public void setBox(Box box) {
        super.setBox(box);
        iterator.setBox(box);

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];

        int totalWV = waveVectors.length;
        
        Qr = new double[totalWV][coordinateDim];
        Qi = new double[totalWV][coordinateDim];
    }

    protected void doStepInternal() {
        super.doStepInternal();
        iterator.reset();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();

        double sqrtT = Math.sqrt(temperature);

        int totalWV = waveVectors.length;
    
   
       	for (int iVector=0; iVector< totalWV; iVector++) {
       
       		for (int j=0; j<coordinateDim; j++) {
       			if (omega[iVector][j] == Double.POSITIVE_INFINITY) continue;
                
       			//generate real and imaginary parts of random normal-mode coordinate Q
       			Qr[iVector][j] = Math.cos(-omega[iVector][j]*currentTime) *sqrtT;
       			Qi[iVector][j] = Math.sin(-omega[iVector][j]*currentTime) *sqrtT;
                
       		}
       	}
       
        
        for (int iCell = 0; iCell<cells.length; iCell++) {

            BasisCell cell = cells[iCell];
            for (int i=0; i<coordinateDim; i++) {
                u[i] = 0;
            }
            //loop over wavevectors and sum contribution of each to the generalized coordinates
            
            /*
             * u up to this point here is the normal mode coordinate
             * 	and then tranformed back to real coordinate in 
             * 	the next part
             */
            if (isOneWV()){
            	int wvNum = getWaveVectorNum();
            	int eValNum = getEValNum();
            	
            	double kR = waveVectors[wvNum].dot(cell.cellPosition);//getLatticePositions()[atomCount]);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                 
                if (isOneEVal()){
                	for (int i=0; i<coordinateDim; i++) {
                		for (int j=0; j<coordinateDim; j++) {
                			u[j] += waveVectorCoefficients[wvNum]*eigenVectors[wvNum][eValNum][j]*
                                  2.0*(Qr[wvNum][eValNum]*coskR - Qi[wvNum][eValNum]*sinkR);
                		}
                	}
                } else {
                	for (int i=0; i<coordinateDim; i++) {
                		for (int j=0; j<coordinateDim; j++) {
                			u[j] += waveVectorCoefficients[wvNum]*eigenVectors[wvNum][i][j]*
                                  2.0*(Qr[wvNum][i]*coskR - Qi[wvNum][i]*sinkR);
                		}
                	}
                }
                
            } else {
            
	            for (int iVector=0; iVector< totalWV; iVector++) {
	            	
	                double kR = waveVectors[iVector].dot(cell.cellPosition);//getLatticePositions()[atomCount]);
	                double coskR = Math.cos(kR);
	                double sinkR = Math.sin(kR);
	                
	                for (int i=0; i<coordinateDim; i++) {
	                    for (int j=0; j<coordinateDim; j++) {
	                        u[j] += waveVectorCoefficients[iVector]*eigenVectors[iVector][i][j]*
	                                  2.0*(Qr[iVector][i]*coskR - Qi[iVector][i]*sinkR);
	                    }
	                }
	            }
            }
            
            double normalization = 1/Math.sqrt(cells.length);
            for (int i=0; i<coordinateDim; i++) {
                u[i] *= normalization;
            }
            coordinateDefinition.setToU(cell.molecules, u);
        }
    }
    
	public void setOneWV(boolean b) {
		oneWV = b;
	}

	public boolean isOneWV() {
		return oneWV;
	}

	public int getWaveVectorNum() {
		return waveVectorNum;
	}

	public void setWaveVectorNum(int waveVectorNum) {
		this.waveVectorNum = waveVectorNum;
	}
	
	public boolean isOneEVal() {
		return oneEVal;
	}

	public void setOneEVal(boolean oneEVal) {
		this.oneEVal = oneEVal;
	}

	public int getEValNum() {
		return eValNum;
	}

	public void setEValNum(int valNum) {
		eValNum = valNum;
	}

	private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final MoleculeIteratorAllMolecules iterator;
    private double[][] omega;
    private double[][][] eigenVectors;
    private Vector[] waveVectors;
    private double[] waveVectorCoefficients;
    protected double[] u;
    protected double[][] Qr;
    protected double[][] Qi;
    protected double lastEnergy;
    protected double temperature;
    protected boolean isRejectable;
    protected double[][] uOld;
    protected boolean oneWV = false, oneEVal = false;
    protected int waveVectorNum, eValNum;

}
