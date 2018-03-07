/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * A Monte Carlo move which compares several normal modes to harmonic normal
 * modes and then explores the phase space defined by the remaining normal
 * modes.
 * 
 * This move has a whitelist of wavevectors that are allowed to be changed, 
 * because they are not harmonic wavevectors.  (This differs from 
 * MCMoveCompareMultipleModes.)
 * 
 * comparedWVs are the wavevectors that are compared (removed).
 * 
 * @author cribbin
 * 
 */
public class MCMoveCompareMultipleWV extends MCMoveBoxStep {

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorLeafAtoms iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld,energyNew, energyEvenLater;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private Vector[] waveVectors;
    private double[] gaussian;
    protected double temperature;
    private double[] rRand, iRand, realT, imagT;
    private double[] waveVectorCoefficients, sqrtWVC;
    private double[][] omega2, oneOverOmega2;
    double[] uNow;
    int changedWV, howManyChangesToHardRodModes;
    int[] comparedWVs, changeableWVs;
    
    public MCMoveCompareMultipleWV(PotentialMaster potentialMaster,
            IRandom random) {
        super(potentialMaster);
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        gaussian = new double[2];
        howManyChangesToHardRodModes = 1;
    }
    
    public boolean doTrial() {
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        rRand = new double[coordinateDim];
        iRand = new double[coordinateDim];
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1 / Math.sqrt(cells.length);
        int numWV = comparedWVs.length;
        
        /*
         * This loop looks at each wavevector, asks if that wavevector is
         * removed, and calculates what happens if it is.
         */
// ZERO OUT A NORMAL MODE.
        
        //Store old positions
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
        }
        
        for (int wvCount = 0; wvCount < numWV; wvCount++) {
            int comparedwv = comparedWVs[wvCount];
            
            // Get normal mode coordinate information
            coordinateDefinition.calcT(waveVectors[comparedwv], realT, imagT);
            
            for (int iCell = 0; iCell < cells.length; iCell++) {
                cell = cells[iCell];
                uNow = coordinateDefinition.calcU(cells[iCell].molecules);
                // rezero deltaU
                for (int j = 0; j < coordinateDim; j++) {
                    deltaU[j] = 0.0;
                }

                // Calculate the contributions to the current position of
                // the zeroed mode, and subtract it from the overall position
                double kR = waveVectors[comparedwv].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                double[][] realCoord = new double[waveVectors.length][coordinateDim];
                double[][] imagCoord = new double[waveVectors.length][coordinateDim];
                for (int iMode = 0; iMode < coordinateDim; iMode++) {// Calculate the current coordinate:
                    realCoord[comparedwv][iMode] = 0.0;
                    imagCoord[comparedwv][iMode] = 0.0;
                    for (int j = 0; j < coordinateDim; j++) {
                        realCoord[comparedwv][iMode] += eigenVectors
                                [comparedwv][iMode][j] * realT[j];
                        imagCoord[comparedwv][iMode] += eigenVectors
                                [comparedwv][iMode][j] * imagT[j];
                    }
                    for (int j = 0; j < coordinateDim; j++) {
                        deltaU[j] -= sqrtWVC[comparedwv] * 
                                eigenVectors[comparedwv][iMode][j] * 
                                (realCoord[comparedwv][iMode] * coskR 
                                - imagCoord[comparedwv][iMode] * sinkR);
                    }
                }
                for (int i = 0; i < coordinateDim; i++) {
                    deltaU[i] *= normalization;
                }

                for (int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }//end of cell loop
        }// end of wvCount loop
        
        energyOld = energyMeter.getDataAsScalar();
        if (Double.isInfinite(energyOld)) {
            int limit = coordinateDefinition.getBox().getLeafList().size();
            for (int k = 0; k < limit; k++) {
                System.out.println(k + " " + 
                        coordinateDefinition.getBox().getLeafList().get(k)
                        .getPosition());
            }
            throw new IllegalStateException(
                    "Overlap after the removal of a mode!");
        }

// MOVE SOME NUMBER OF RANDOM HARD ROD POTENTIAL MODES, AND MEASURE energyNew
        // equivalent to MCMoveChangeMode for several modes.
        for(int wvCount = 0; wvCount < howManyChangesToHardRodModes; wvCount++){
            
            // Select the wave vector whose eigenvectors will be changed.
            changedWV = changeableWVs[random.nextInt(changeableWVs.length)];
            
            // calculate the new positions of the atoms.
            // loop over cells
            double[] delta = new double[coordinateDim*2];
            for ( int i = 0; i < coordinateDim*2; i++) {
                delta[i] = (2*random.nextDouble()-1) * stepSize;
            }
            for (int iCell = 0; iCell < cells.length; iCell++) {
                uNow = coordinateDefinition.calcU(cells[iCell].molecules);
                cell = cells[iCell];
                // rezero deltaU
                for (int j = 0; j < coordinateDim; j++) {
                    deltaU[j] = 0.0;
                }
                // loop over the wavevectors, and sum contribution of each to
                // the generalized coordinates. Change the selected wavevectors
                // eigenvectors at the same time!
                double kR = waveVectors[changedWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int iMode = 0; iMode < coordinateDim; iMode++) {
                    if( !(Double.isInfinite(omega2[changedWV][iMode]))) {
                        for (int j = 0; j < coordinateDim; j++) {
                            deltaU[j] += sqrtWVC[changedWV]*eigenVectors[changedWV][iMode][j]*
                                (delta[iMode]*coskR - delta[iMode+coordinateDim]*sinkR);
                        }
                    }
                }
                 for(int i = 0; i < coordinateDim; i++){
                     deltaU[i] *= normalization;
                 }
                for (int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }
        }//end wvCount loop
        energyNew = energyMeter.getDataAsScalar();
        
//        System.out.println("after hardrod move: " + energyNew);
//        for(int i = 0; i < list.getAtomCount(); i++){
//            System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//        }
        
// MOVE EACH NORMAL MODE THAT WAS ZEROED OUT.
        // set up the gaussian values
        double sqrtT = Math.sqrt(temperature);
        
        //This should loop over the wave vectors that we are comparing.
        for(int wvCount = 0; wvCount < numWV; wvCount++){
            int comparedwv = comparedWVs[wvCount];
            
            for (int j = 0; j < coordinateDim; j++) {
                if (oneOverOmega2[comparedwv][j] == 0) {continue;}
                // generate real and imaginary parts of random normal-emode
                // coordinate Q
                double realGauss = random.nextGaussian() * sqrtT;
                double imagGauss = random.nextGaussian() * sqrtT;
                
                // XXX we know that if c(k) = 0.5, one of the gaussians will be
                // ignored, but it's hard to know which. So long as we don't put
                // an atom at the  origin (which is true for 1D if c(k)=0.5), 
                // it's the real part that will be ignored.
                if (waveVectorCoefficients[comparedwv] == 0.5) {imagGauss = 0;}
                rRand[j] = realGauss * oneOverOmega2[comparedwv][j];
                iRand[j] = imagGauss * oneOverOmega2[comparedwv][j];
                gaussian[0] = realGauss;
                gaussian[1] = imagGauss;
            }
    
            // calculate the new positions of the atoms.
            for (int iCell = 0; iCell < cells.length; iCell++) {
                uNow = coordinateDefinition.calcU(cells[iCell].molecules);
                cell = cells[iCell];
                // rezero deltaU
                for (int j = 0; j < coordinateDim; j++) {
                    deltaU[j] = 0.0;
                }
                // Calculate the change in position due to the substitution of a
                // Gaussian.
                double kR = waveVectors[comparedwv].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int i = 0; i < coordinateDim; i++) {
                    for (int j = 0; j < coordinateDim; j++) {
                        deltaU[j] += sqrtWVC[comparedwv] * eigenVectors[comparedwv][i][j] * 
                                 (rRand[i] * coskR - iRand[i] * sinkR);
                    }
                }
                for (int i = 0; i < coordinateDim; i++) {
                    deltaU[i] *= normalization;
                }
                for (int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }
                
        } // end wvCount loop
        
//        System.out.println("At end of move: " + energyMeter.getDataAsScalar());
//        for(int i = 0; i < list.getAtomCount(); i++){
//            System.out.println(((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(i)).getPosition());
//        }
        return true;
    }

    public double getChi(double temperature) {
        return Math.exp(-(energyNew - energyOld) / temperature);
    }

    public void acceptNotify() {
//      System.out.println("accept MCMoveCompareMultipleWV");
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//      System.out.println("reject MCMoveCompareMultipleWV");
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell < cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);
        energyMeter.setBox(newBox);
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public void setCoordinateDefinition(CoordinateDefinition newCD) {
        coordinateDefinition = newCD;
        deltaU = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;
        realT = new double[coordinateDefinition.getCoordinateDim()];
        imagT = new double[coordinateDefinition.getCoordinateDim()];
    }

    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    /**
     * Set the wave vectors accessible to the move.
     * 
     * @param wv
     */
    public void setWaveVectors(Vector[] wv){
        waveVectors = new Vector[wv.length];
        waveVectors = wv;
    }
    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        waveVectorCoefficients = newWaveVectorCoefficients;
        
        sqrtWVC = new double[newWaveVectorCoefficients.length];
        for (int i =0; i < newWaveVectorCoefficients.length; i++){
            sqrtWVC[i] = Math.sqrt(2*newWaveVectorCoefficients[i]);
        }
    }

    /**
     * Informs the move of the eigenvectors
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }

    public void setOmegaSquared(double[][] o2, double[] coeff) {
        this.omega2 = o2;
        oneOverOmega2 = new double[o2.length][o2[0].length];
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = Math.sqrt(1.0/(omega2[i][j]));
            }
        }
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    
    }

    public double[] getGaussian() {
        return gaussian;
    }

    /**
     * Set the wavevectors that are removed.
     * @param wv the wavevectors that are removed.
     */
    public void setComparedWV(int[] wv){
        comparedWVs = wv;
    
        int dink = 0;
        for(int i = 0; i < wv.length; i++){
            if(wv[i] < wv[dink]){
                dink = i;
            }
        }
    }
    
    /**
     * Set the wavevectors that can be changed by the move.  Compared wavevectors
     * should not be on this list.
     * @param wv the wavevectors that can be changed by the move.
     */
       
    public void setChangeableWVs(int[] wv){
        changeableWVs = wv;
        
    }
    
    private void printLocations(){
        IAtomList list = box.getLeafList();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int ats = box.getLeafList().size();
        
        if(box.getBoundary().getEdgeVector(0).getD() == 1){
            for(int i = 0; i < ats; i++){
                System.out.println(i + "  " + list.get(i).getPosition().getX(0));
            }
        }
        
        if(box.getBoundary().getEdgeVector(0).getD() == 3){
            for(int i = 0; i < ats; i++){
                System.out.println("Atom " + i);
                for(int j = 0; j < 3; j++){
                    System.out.println(j + " " + list.get(i).getPosition().getX(j));
                }
            }
        }
    }
}
