/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.api.*;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;

/**
 * A Monte Carlo move which selects a wave vector, and changes the normal modes
 * associated with that wave vector.
 * 
 * Uses a whitelist of wavevectors that can be changed by the doTrial() method.
 * 
 * @author cribbin
 *
 */
public class MCMoveChangeMultipleWV extends MCMoveBoxStep{

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorLeafAtoms iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld, energyNew;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private double[][] omega2, oneOverOmega2;
    private Vector[] waveVectors;
    private double[] waveVectorCoefficients, sqrtWVC;
    int changedWV;
    int[] changeableWV;  //all wvs from this are changed.
    
    
    public MCMoveChangeMultipleWV(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        
        changeableWV = new int[0];
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        deltaU = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }
    
    public void addChangeableWV(int[] wv){
        
        int[] newWV = new int[wv.length + changeableWV.length];
        
        for (int i = 0; i < wv.length; i++){
            for (int j = 0; j < changeableWV.length; j++){
                if(wv[i] == changeableWV[j]){
                    throw new IllegalArgumentException("WV used twice in " +
                            "MCMoveChangeMultipleWV.addChangeableWV");
                }
            }
            newWV[i] = wv[i];
        }
        
        for (int i = wv.length; i < newWV.length; i++){
            newWV[i] = changeableWV[i-wv.length];
        }
        
        changeableWV = new int[newWV.length];
        changeableWV = newWV;
        
    }
    
    public void addChangeableWV(String string){
      changeableWV = new int[waveVectors.length];
        for(int i = 0; i < changeableWV.length; i++){
            changeableWV[i] = i;
        }
    }
    
    public void setOmegaSquared(double[][] o2){
        omega2 = o2;
        oneOverOmega2 = new double[o2.length][o2[0].length];
        for (int i=0; i<oneOverOmega2.length; i++) {
            for (int j=0; j<oneOverOmega2[i].length; j++) {
                oneOverOmega2[i][j] = Math.sqrt(1.0/(omega2[i][j]));
            }
        }
    }
    /**
     * Set the wave vectors that can be changed by the move.
     * 
     * @param wv
     */
    public void setWaveVectors(Vector[] wv){
        waveVectors = new Vector[wv.length];
        waveVectors = wv;
    }
    public void setWaveVectorCoefficients(double[] coeff){
        waveVectorCoefficients = coeff;
        sqrtWVC = new double[coeff.length];
        for (int i =0; i < coeff.length; i++){
            sqrtWVC[i] = Math.sqrt(2*coeff[i]);
        }
    }
    /**
     * Informs the move of the eigenvectors for the selected wave vector.  The
     * actual eigenvectors used will be those specified via setModes
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
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
        
        energyOld = energyMeter.getDataAsScalar();
        
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        
        // assume that the first cell is the same as every other cell.
        //Also, this sets up the variable "cell".
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        
        // Select the wave vector whose eigenvectors will be changed.
        changedWV = changeableWV[random.nextInt(changeableWV.length)];
        
        //calculate the new positions of the atoms.
        //loop over cells
        double[] delta = new double[coordinateDim*2];
        
        for (int i = 0; i < coordinateDim*2; i++){
            delta[i] = (2*random.nextDouble()-1) * stepSize;
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
            
            //N.B. i is the mode, j is the coordinate (atom+direction (3D))
            for(int iMode = 0; iMode < coordinateDim; iMode++){
            
                if( !(Double.isInfinite(omega2[changedWV][iMode])) ){
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] += sqrtWVC[changedWV]*eigenVectors[changedWV][iMode][j]*
                            (delta[iMode]*coskR - delta[iMode+coordinateDim]*sinkR);
                    }
                }
            }
            double normalization = 1/Math.sqrt(cells.length);
            for(int i = 0; i < coordinateDim; i++){
                deltaU[i] *= normalization;
            }
            
            for(int i = 0; i < coordinateDim; i++) {
                uNow[i] += deltaU[i];
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            
        }//end of cell loop
        
        energyNew = energyMeter.getDataAsScalar();
        
//        System.out.println("energy new " + energyNew);
//        printLocations();

        return true;
    }
    
    public double getA() {
        return 1;
    }

    public double getB() {
        return -(energyNew - energyOld);
    }
    
    public void acceptNotify() {
//        System.out.println("accept MCMoveChangeMultipleWV");
//        System.out.println("nrg dif" + (energyOld - energyNew));
//        printLocations();
        
        iterator.reset();
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//        System.out.println("reject MCMoveChangeMultipleWV");
//        System.out.println("nrg dif" + (energyOld - energyNew));
//        printLocations();
        
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }

    private void printLocations(){
        IAtomList list = box.getLeafList();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int ats = box.getLeafList().getAtomCount();
        
        if(box.getBoundary().getEdgeVector(0).getD() == 1){
            for(int i = 0; i < ats; i++){
                System.out.println(i + "  " + list.getAtom(i).getPosition().getX(0));
            }
        }
        
        if(box.getBoundary().getEdgeVector(0).getD() == 3){
            for(int i = 0; i < ats; i++){
                System.out.println("Atom " + i);
                for(int j = 0; j < 3; j++){
                    System.out.println(j + " " + list.getAtom(i).getPosition().getX(j));
                }
            }
        }
    }

}
