/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.IPotential;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.units.Null;


/**
 * This meter assumes the use of only MCMoveCompareMode in the simulation.
 * Assumes 1D system - otherwise, choose a mode and eliminate i loops.
 * 
 * @author cribbin
 *
 */
public class MeterCompareWVShortcut extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    
    private double eigenVectors[][][];
    private Vector[] waveVectors;
    int comparedWV;
    protected double temperature;
    private double[] waveVectorCoefficients, sqrtWVC;
    private double wvc;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double[][] omegaSquared, oneOverOmega2; //spring constants
    private double energyNM, energyOP;  //energyNM is energy of normal modes
                                //energyOP is the energy of the Gaussian modes
    
    private static final long serialVersionUID = 1L;
    MCMoveCompareMultipleWV mcmove;
    
    public MeterCompareWVShortcut(PotentialMaster potentialMaster,
            MCMoveCompareMultipleWV mcmove){
        super("meterComapreMode", Null.DIMENSION);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        
        this.mcmove = mcmove; 
    }
    public double getDataAsScalar() {
        double[] gaussian = mcmove.getGaussian();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        energyNM = 0.0;
        energyOP = 0.0;
        
        //get normal mode coordinate of "last" wavevector
        coordinateDefinition.calcT(waveVectors[comparedWV], realT, imagT);
        double realCoord = 0.0, imagCoord = 0.0;
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            for(int j = 0; j < coordinateDim; j++){
                realCoord += eigenVectors[comparedWV][i][j] * realT[j];
                imagCoord += eigenVectors[comparedWV][i][j] * imagT[j];
            }
        }
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store original positions
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            for(int j = 0; j < coordinateDim; j++){
                deltaU[j] = 0.0;
            }
            
            //Calculate the contributions to the current position of the 
            //zeroed mode, and subtract it from the overall position.
            double kR = waveVectors[comparedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                //Calculate the current coordinates
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= sqrtWVC[comparedWV]*eigenVectors[comparedWV][i][j] *
                        (realCoord*coskR - imagCoord*sinkR);
                }
            }
            
            for(int i = 0; i < coordinateDim; i++) {  //Loop would go away
                uNow[i] += deltaU[i];
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
        }
        energyNM = meterPE.getDataAsScalar();
        
        //Calculate the energy due to the Gaussian modes.
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            if(Double.isInfinite(omegaSquared[comparedWV][i])){
                continue;
            }
          double normalCoord = gaussian[0]*gaussian[0] + gaussian[1]*gaussian[1];
            energyOP += wvc * normalCoord * omegaSquared[comparedWV][i];
        }
        
        // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
        return energyNM + energyOP;
    }

    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(Vector[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setComparedWV(int comparedWV) {
        this.comparedWV = comparedWV;
        wvc = waveVectorCoefficients[comparedWV];
    }
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    public void setWaveVectorCoefficients(double[] waveVectorCoefficients) {
        this.waveVectorCoefficients = waveVectorCoefficients;
        sqrtWVC = new double[waveVectorCoefficients.length];
        for (int i =0; i < waveVectorCoefficients.length; i++){
            sqrtWVC[i] = Math.sqrt(2*waveVectorCoefficients[i]);
        }
    }
    public void setCoordinateDefinition(CoordinateDefinition cd){
        coordinateDefinition = cd;
        coordinateDim = coordinateDefinition.getCoordinateDim();
    }
    public void setSpringConstants(double[][] sc){
        omegaSquared = sc;
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = Math.sqrt(1.0/(sc[i][j]));
            }
        }
    }
    public void setOmegaSquared(double[][] sc){
        omegaSquared = sc;
        oneOverOmega2 = new double[sc.length][sc[0].length];
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = Math.sqrt(1.0/(sc[i][j]));
            }
        }
    }
    
    

}
