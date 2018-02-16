/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.dimensions.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * 
 * @author cribbin
 *
 */
public class MeterWidomMode extends DataSourceScalar {

    public int nInsert, affectedWV, counter;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition coordinateDefinition;
    private int coordinateDim;
    private double eigenVectors[][][];
    private Vector[] waveVectors;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared, oneOverOmega2;
    protected double temperature;
    private double[] uNow, deltaU;
    private double[] sqrtWVC;
    
    
    public MeterWidomMode(String string, PotentialMaster
            potentialMaster, CoordinateDefinition cd, Box box, int awv){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster, box);
        affectedWV = awv;
    }
    
    public double getDataAsScalar() {
        
//      IAtomList leaflist = coordinateDefinition.getBox().getLeafList();
//      System.out.println("start:");
//      for(int i = 0; i < 32; i++){
//          System.out.println(i + "  " + ((Atom)leaflist.getAtom(i)).getPosition().x(0));
//      }
        
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        
        //get normal mode coordinate of affected waveVector
        coordinateDefinition.calcT(waveVectors[affectedWV], realT, imagT);
        
        double[] realCoord = new double[coordinateDim];
        double[] imagCoord = new double[coordinateDim];
        for(int j = 0; j < coordinateDim; j++){
            realCoord[j] = 0.0;
            imagCoord[j] = 0.0;
        }
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            if(Double.isInfinite(omegaSquared[affectedWV][i])){
                continue;
            }// end if
            for(int j = 0; j < coordinateDim; j++){
                realCoord[i] += eigenVectors[affectedWV][i][j] * realT[j];
                imagCoord[i] += eigenVectors[affectedWV][i][j] * imagT[j];
            }
        }//end loops which sum realCoord and imagCoord
         
        //remove affected wavevector's effect on position
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store original positions
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            for(int j = 0; j < coordinateDim; j++){
                deltaU[j] = 0.0;
            }//end zero-the-deltaUs-loop
            
            //Calculate the contributions to the current position of the 
            //zeroed mode, and subtract it from the overall position.
            double kR = waveVectors[affectedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinates.
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= sqrtWVC[i]*eigenVectors[affectedWV][i][j] *
                        (realCoord[i]*coskR - imagCoord[i]*sinkR);
                }
            }

            for(int i = 0; i < coordinateDim; i++){
                deltaU[i] *= normalization;
            }
            for(int i = 0; i < coordinateDim; i++) {
                uNow[i] += deltaU[i];
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
        }//end of cell loop
        
        //Calculate the energy without that mode in.
        double energy = meterPE.getDataAsScalar();
        
        // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
        if(Double.isInfinite(energy)) {
            return 0;
        } else {
            return 1;
        }
    }

    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(Vector[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setAffectedWV(int awv) {
        this.affectedWV = awv;
    }
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    public void setWaveVectorCoefficients(double[] waveVectorCoefficients) {
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
        oneOverOmega2 = new double[sc.length][sc[0].length];
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
