/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.simulation.Simulation;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.NormalModes;
import etomica.space.Space;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume a mode & rod 
 * are removed
 * 
 * Should only be used when one system is exactly twice the size of the other
 * system.
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageSubtractDoubleSize extends MeterDifferentImageSubtract {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MeterDifferentImageSubtractDoubleSize";
    
    public MeterDifferentImageSubtractDoubleSize(Simulation sim, Space space,
                                                 CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells,
                                                 NormalModes otherNM){
        this(sim, space, simCD, simNM, otherCD, potentialMaster, 
                otherNCells, otherNM, "file");
    }
    public MeterDifferentImageSubtractDoubleSize(Simulation sim, Space space,
                                                 CoordinateDefinition simCD, NormalModes simNM, CoordinateDefinition
            otherCD, PotentialMasterList potentialMaster, int[] otherNCells,
                                                 NormalModes otherNM, String otherFilename){
        super(sim, space, simCD, simNM, otherCD, potentialMaster, otherNCells,
                otherNM, otherFilename);
    }
    
    public double getDataAsScalar() {
        double normalization = 1/Math.sqrt(2.0);  //DS
        
        //Calculate normal mode coordinates of simulation system.
        // First index is wavevector, second index is mode.
        double[] realCoord = new double[simCDim];
        double[] imagCoord = new double[simCDim];
        
        simCDef.calcT(simWaveVectors[1], simRealT, simImagT);
        for (int iMode = 0; iMode < simCDim; iMode++){
            realCoord[iMode] = 0.0;
            imagCoord[iMode] = 0.0;
            for (int j = 0; j < simCDim; j++){
                realCoord[iMode] += simEigenVectors[1][iMode][j] * simRealT[j];
                imagCoord[iMode] += simEigenVectors[1][iMode][j] * simImagT[j];
            }
            if(simWVCoeff[1] == 1.0){
                realCoord[iMode] *= Math.sqrt(2);
                imagCoord[iMode] *= Math.sqrt(2);
            }
        }
        
        
        //CALCULATE A HARMONIC ENERGY FOR THE SUBTRACTED NORMAL MODES.
        //Calculate harmonic, and zero out the subtracted normal modes.
        double harmonic = 0.0;
        if(simWVCoeff[1] == 1.0){
            for (int iMode = 0; iMode < simCDim; iMode++){
                if(!(sqrtSimOmega2[1][iMode] == Double.POSITIVE_INFINITY)){
                    harmonic += 0.5 * realCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * realCoord[iMode] * sqrtSimOmega2[1][iMode];
                    harmonic += 0.5 * imagCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * imagCoord[iMode] * sqrtSimOmega2[1][iMode];
                }
            }
        } else {
            for (int iMode = 0; iMode < simCDim; iMode++){
                if(!(sqrtSimOmega2[1][iMode] == Double.POSITIVE_INFINITY)){
                    harmonic += 0.5 * realCoord[iMode] * sqrtSimOmega2[1][iMode]
                        * realCoord[iMode] * sqrtSimOmega2[1][iMode];
                }
            }
        }
        
        
        //Calculate the positions for the meter's system, and set them.
        double[] temp = simCDef.calcU(simCDef.getBasisCells()[0].molecules);
        double[] simZero = new double[temp.length];
        double[] simOne = new double[temp.length];
        for(int i = 0; i < temp.length; i++){
            simZero[i] = temp[i];
        }
        temp = simCDef.calcU(simCDef.getBasisCells()[1].molecules);
        for(int i = 0; i < temp.length; i++){
            simOne[i] = temp[i];
        }
        double[] newU = new double[simZero.length];
        for (int i = 0; i < simZero.length; i++){
            newU[i] = (simZero[i]+ simOne[i]);
            newU[i] *= normalization;
        }
        cDef.setToU(cDef.getBasisCells()[0].molecules, newU);
        
        return meterPE.getDataAsScalar() + harmonic;
    }

}
