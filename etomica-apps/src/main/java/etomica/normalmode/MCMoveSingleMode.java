/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.space.Vector;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition.BasisCell;

/**
 * A Monte Carlo move which selects each individual modes to move
 * the atoms accordingly
 * 
 * @author Tai Boon Tan
 *
 */
public class MCMoveSingleMode extends MCMoveBoxStep{

    public MCMoveSingleMode(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        deltaU = new double[coordinateDefinition.getCoordinateDim()];
        delta1 = new double[coordinateDefinition.getCoordinateDim()];
        delta2 = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }
    
    /**
     * The harmonic wavevector and all wavevectors with higher numbers are not
     * able to be changed by this MCMove.
     */
    public void setHarmonicWV(int hwv){
        harmonicWV = hwv;
    }
    
    public void setTemperature(double temp){
    	temperature = temp;
    }
    
    /**
     * Set the wave vectors used by the move.
     * 
     * @param wv
     */
    public void setWaveVectors(Vector[] wv){
        waveVectors = new Vector[wv.length];
        waveVectors = wv;
    }
    public void setWaveVectorCoefficients(double[] coeff){
        waveVectorCoefficients = coeff;
    }
    /**
     * Informs the move of the eigenvectors for the selected wave vector.  The
     * actual eigenvectors used will be those specified via setModes
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] o2, double[] coeff){
    	omega2 = o2;
    	/*
    	scaledStepSize = new double[omega2.length][omega2[0].length];
    	for (int i=0; i<scaledStepSize.length; i++) {
            for (int j=0; j<scaledStepSize[i].length; j++) {
            	scaledStepSize[i][j] = Math.sqrt(1.0/(2.0*omega2[i][j]*coeff[i]));
            }
        }
        */
    }
    
    public void setBox(Box newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);
        energyMeter.setBox(newBox);
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    
    public boolean doTrial() {
//        System.out.println("initial positions:");
//        printLocations();
        
        energyOld = energyMeter.getDataAsScalar();
//        System.out.println("energyOld " +energyOld);
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        
        BasisCell cell;
//        double[] calcedU = coordinateDefinition.calcU(cell.molecules);
        uOld = new double[cells.length][coordinateDim];
        
        // Select the wave vector whose eigenvectors will be changed.
        //The zero wavevector is center of mass motion, and is rejected as a 
        //possibility.
        changedWV = random.nextInt(harmonicWV);
        
        //calculate the new positions of the atoms.
        //loop over cells
        for (int i=0; i<coordinateDim; i++){
        	delta1[i] = (2*random.nextDouble()-1) *stepSize;
        	delta2[i] = (2*random.nextDouble()-1) *stepSize;
        }
        
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            double[] uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            for(int i = 0; i< coordinateDim; i++){
                  deltaU[i] = 0;
            }
            
            //loop over the wavevectors, and sum contribution of each to the
            //generalized coordinates.  Change the selected wavevectors eigen-
            //vectors at the same time!
            double kR = waveVectors[changedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                if( !(Double.isInfinite(omega2[changedWV][i])) ){
                    //if(changedWV == 0) { System.out.println("why am i here?");}
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] += waveVectorCoefficients[changedWV] * 
                            eigenVectors[changedWV][i][j]*2.0*(delta1[i]*coskR - delta2[i]*sinkR);
                    }
                }
            }
            double normalization = 1/Math.sqrt(cells.length);
            for(int i = 0; i < coordinateDim; i++){
                deltaU[i] *= normalization;
            }
            
            for(int i = 0; i < coordinateDim; i++) {
                uNow[i] += deltaU[i];
                //System.out.println("uNow["+i+"]: "+ uNow[i]);
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            
        }
        
        energyNew = energyMeter.getDataAsScalar();
//        System.out.println("energyNew " + energyNew);
        return true;
    }
    
    public double getA() {
        return 1;
    }

    public double getB() {
        return -(energyNew - energyOld);
    }
    
    public void acceptNotify() {
//        System.out.println("accept MCMoveChangeSingleMode");
//        printLocations();
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//        System.out.println("reject MCMoveChangeSingleMode");
//        printLocations();
        
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
//        System.out.println("After moved back");
//        printLocations();
    }

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorLeafAtoms iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld, energyNew /*, latticeEnergy*/;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private double[][] omega2, scaledStepSize;
    private Vector[] waveVectors;
    private double[] waveVectorCoefficients;
    int changedWV, harmonicWV;  //all wvs from the harmonic wv and up are not changed.
    protected double temperature;
    protected double[] delta1, delta2;
    
    
}
