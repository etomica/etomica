/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveHarmonic extends MCMoveBox {

    public MCMoveHarmonic(IRandom random) {
        super(null, new MCMoveTracker());
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
    }

    public void setRejectable(boolean newIsRejectable) {
        isRejectable = newIsRejectable;
    }
    
    public boolean isRejectable() {
        return isRejectable;
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setOmegaSquared(double[][] omega2) {
        oneOverOmega2 = new double[omega2.length][omega2[0].length];
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = 1.0/Math.sqrt(omega2[i][j]);
            }
        }
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        sqrtWVC = new double[newWaveVectorCoefficients.length];
        for (int i =0; i < newWaveVectorCoefficients.length; i++){
            sqrtWVC[i] = Math.sqrt(2*newWaveVectorCoefficients[i]);
        }
    }
    
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setBox(Box newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];

        rRand = new double[waveVectors.length][coordinateDim];
        iRand = new double[waveVectors.length][coordinateDim];
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public boolean doTrial() {
        iterator.reset();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();

        lastEnergy = 0;
        double sqrtT = Math.sqrt(temperature);

        for (int iVector=0; iVector<waveVectors.length; iVector++) {
            
            if (isSelectMode){
               for (int j=0; j<modeNum.length; j++) {
                   if (oneOverOmega2[iVector][modeNum[j]] == 0) continue;
                   //generate real and imaginary parts of random normal-mode coordinate Q
                   double realGauss = random.nextGaussian() * sqrtT;
                   double imaginaryGauss = random.nextGaussian() * sqrtT;
                   rRand[iVector][modeNum[j]] = realGauss * oneOverOmega2[iVector][modeNum[j]];
                   iRand[iVector][modeNum[j]] = imaginaryGauss * oneOverOmega2[iVector][modeNum[j]];
                   //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
                   // it's hard to know which.  So long as we don't put an atom at the origin
                   // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
                   if (sqrtWVC[iVector] == 1.0) imaginaryGauss = 0;
                   lastEnergy += 0.5 * (realGauss*realGauss + imaginaryGauss*imaginaryGauss);
               }
            } else {
                for (int j=0; j<coordinateDim; j++) {
                    if (oneOverOmega2[iVector][j] == 0) continue;
                    //generate real and imaginary parts of random normal-mode coordinate Q
                    double realGauss = random.nextGaussian() * sqrtT;
                    double imaginaryGauss = random.nextGaussian() * sqrtT;
                    rRand[iVector][j] = realGauss * oneOverOmega2[iVector][j];
                    iRand[iVector][j] = imaginaryGauss * oneOverOmega2[iVector][j];
                    //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
                    // it's hard to know which.  So long as we don't put an atom at the origin
                    // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
                    if (sqrtWVC[iVector] == 1.0) imaginaryGauss = 0;
                    lastEnergy += 0.5 * (realGauss*realGauss + imaginaryGauss*imaginaryGauss);
                }
            }
        }
        
        if (isRejectable) {
            if (uOld == null || uOld.length != cells.length) {
                uOld = new double[cells.length][coordinateDim];
            }
        }
        
        double normalization = 1/Math.sqrt(cells.length);
        for (int iCell = 0; iCell<cells.length; iCell++) {
            if (isRejectable) {
                double[] uNow = coordinateDefinition.calcU(cells[iCell].molecules);
                System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            }
            BasisCell cell = cells[iCell];
            for (int i=0; i<coordinateDim; i++) {
                u[i] = 0;
            }
            //loop over wavevectors and sum contribution of each to the generalized coordinates
            for (int iVector=0; iVector<waveVectors.length; iVector++) {
                double kR = waveVectors[iVector].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                
                if(isSelectMode){
                    for (int i=0; i<modeNum.length; i++) {
                        for (int j=0; j<coordinateDim; j++) {
                            u[j] += sqrtWVC[iVector]*eigenVectors[iVector][modeNum[i]][j]*
                                    (rRand[iVector][modeNum[i]]*coskR - iRand[iVector][modeNum[i]]*sinkR);
                        }
                    }
                } else {
                    for (int i=0; i<coordinateDim; i++) {
                        for (int j=0; j<coordinateDim; j++) {
                            u[j] += sqrtWVC[iVector]*eigenVectors[iVector][i][j]*
                                    (rRand[iVector][i]*coskR - iRand[iVector][i]*sinkR);
                        }
                    }
                }
            }
            for (int i=0; i<coordinateDim; i++) {
                u[i] *= normalization;
            }
            coordinateDefinition.setToU(cell.molecules, u);
        }
        return true;
    }

    public double getChi(double temperature) {
        // return 1 to guarantee success
        return 1;
    }

    /**
     * Returns the harmonic energy of the configuration based on the last
     * harmonic move made by this MC Move.
     */
    public double getLastTotalEnergy() {
        return lastEnergy;
    }
    
    public void acceptNotify() {
    }

    public double energyChange() {
        return 0;
    }

    public void rejectNotify() {
        if (!isRejectable) {
            throw new RuntimeException("I didn't keep track of the old positions.  You have to call setRejectable.");
        }
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }

    
    public int[] getModeNum() {
        return modeNum;
    }

    public void setModeNum(int[] modeNum) {
        isSelectMode = true;
        this.modeNum = modeNum;
    }
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected final AtomIteratorLeafAtoms iterator;
    private double[][] oneOverOmega2;
    private double[][][] eigenVectors;
    private Vector[] waveVectors;
    private double[] sqrtWVC;
    protected double[] u;
    protected double[][] rRand;
    protected double[][] iRand;
    protected final IRandom random;
    protected double lastEnergy;
    protected double temperature;
    protected boolean isRejectable;
    protected double[][] uOld;
    protected int[] modeNum;
    protected boolean isSelectMode = false;

}
