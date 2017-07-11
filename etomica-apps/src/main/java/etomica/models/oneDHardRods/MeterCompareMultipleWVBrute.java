/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.IPotential;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 * Meter which measures a hard rod energy plus a harmonic energy of a system
 * The harmonic energy is determined by removing the effects of some 
 * normal mode coordinates, and applying 1/2 omega^2 nmc^2 to them.
 * The hard rod energy is determined by the remaining wavevectors.
 * @author cribbin
 *
 */
public class MeterCompareMultipleWVBrute extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    private double eigenVectors[][][];
    private Vector[] waveVectors;
    int[] comparedWVs;
    protected double temperature;
    private double[] waveVectorCoefficients, sqrtWVC;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared, oneOverOmega2;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double energyHardRod, energyHarmonic;
    private static final long serialVersionUID = 1L;
    public boolean isOnlyHardRod;
    
    public MeterCompareMultipleWVBrute(PotentialMaster potentialMaster,
            CoordinateDefinition cd, Box box){
        this("meterCompareMultipleModes", potentialMaster, cd, box);
    }
    
    public MeterCompareMultipleWVBrute(String string, PotentialMaster
            potentialMaster, CoordinateDefinition cd, Box box){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
    }
    
    public double getDataAsScalar(){
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        int numWV = comparedWVs.length;
        double[][] realCoord = new double[numWV][coordinateDim];
        double[][] imagCoord = new double[numWV][coordinateDim];
        energyHardRod = 0.0;
        energyHarmonic = 0.0;
        
//        System.out.println("At the start of the meter:");
//        System.out.println("Energy: " + meterPE.getDataAsScalar());
//        IAtomList list = coordinateDefinition.getBox().getLeafList();
//        for(int i = 0; i < list.getAtomCount(); i++){
//            System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//        }
        
        //Get the normal mode coordinates of the compared waveVectors, and
        // store them in realCoord and imagCoord for further use.
        for(int wvcount = 0; wvcount < numWV; wvcount++){
            coordinateDefinition.calcT(waveVectors[comparedWVs[wvcount]], realT, imagT);
            for(int iMode = 0; iMode < coordinateDim; iMode++){  //Loop would go away
                realCoord[wvcount][iMode] = 0.0;
                imagCoord[wvcount][iMode] = 0.0;
                for(int j = 0; j < coordinateDim; j++){
                    realCoord[wvcount][iMode] += eigenVectors[comparedWVs[wvcount]][iMode][j] * realT[j];
                    imagCoord[wvcount][iMode] += eigenVectors[comparedWVs[wvcount]][iMode][j] * imagT[j];
                }
            }
        }
        
        //Remove the effects of the compared modes.
        for(int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            
            //store the original positions of the cell.
            uNow = coordinateDefinition.calcU(cell.molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            
            for(int wvcount = 0; wvcount < numWV; wvcount++){
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] = 0.0;
                }
                
                //Calculate the contributions to the current position of the 
                //zeroed mode, and subtract it from the overall position.
                double kR = waveVectors[comparedWVs[wvcount]].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int iMode = 0; iMode < coordinateDim; iMode++){  //Loop would go away
                    //Calculate the current coordinates.
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] -= sqrtWVC[comparedWVs[wvcount]] * 
                                eigenVectors[comparedWVs[wvcount]][iMode][j] * 
                                (realCoord[wvcount][iMode] * coskR - 
                                imagCoord[wvcount][iMode] * sinkR);
                    }
                }

                for(int i = 0; i < coordinateDim; i++){
                    deltaU[i] *= normalization;
                }
                
                for(int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cell.molecules, uNow);
            }//end of wvcount loop
        }//end of cell loop
        energyHardRod = meterPE.getDataAsScalar();
        
//        System.out.println("After the modes have been removed:");
//        System.out.println("Energy: " + meterPE.getDataAsScalar());
//        list = coordinateDefinition.getBox().getLeafList();
//        for(int i = 0; i < list.getAtomCount(); i++){
//            System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//        }
        
//        if(((Double)energyHardRod).isInfinite() && !isOnlyHardRod) {
//            IAtomList crashlist = coordinateDefinition.getBox().getLeafList();
//            for(int i = 0; i < crashlist.getAtomCount(); i++){
//                System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//            }
//            System.out.println();
//            throw new IllegalArgumentException("Bailing - bad overlap");
//        }
        
        
        //Calculate the energy due to the compared modes
        for(int wvcount = 0; wvcount < numWV; wvcount++){
            for(int iMode = 0; iMode < coordinateDim; iMode++){
                if(Double.isInfinite(omegaSquared[comparedWVs[wvcount]][iMode])){
                    continue;
                }
                double normalCoord = realCoord[wvcount][iMode]*realCoord[wvcount][iMode] + 
                    imagCoord[wvcount][iMode] * imagCoord[wvcount][iMode];
                energyHarmonic += waveVectorCoefficients[comparedWVs[wvcount]] * normalCoord 
                    * omegaSquared[comparedWVs[wvcount]][iMode];
            }
        }
        
        
     // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }

//        System.out.println("At the end of the meter:");
//        System.out.println("Energy: " + meterPE.getDataAsScalar());
//        list = coordinateDefinition.getBox().getLeafList();
//        for(int i = 0; i < list.getAtomCount(); i++){
//            System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//        }
        
        return energyHardRod + energyHarmonic;
    }
    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(Vector[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setComparedWV(int[] cwvs) {
        this.comparedWVs = cwvs;
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

//    public boolean isOnlyHardRod() {
//        return isOnlyHardRod;
//    }
//
//    public void setIsOnlyHardRod(boolean isA) {
//        this.isOnlyHardRod = isA;
//    }
}
