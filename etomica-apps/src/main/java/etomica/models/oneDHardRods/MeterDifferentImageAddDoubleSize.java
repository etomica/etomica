/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.simulation.Simulation;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.normalmode.NormalModes;
import etomica.space.Space;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume an extra mode & rod 
 * are added
 * 
 * Should only be used when one system is exactly twice the size of the other
 * system.
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageAddDoubleSize extends MeterDifferentImageAdd {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MeterDifferentImageAddDoubleSize";
    
    public MeterDifferentImageAddDoubleSize(Simulation sim, Space space,
                                            double temp, CoordinateDefinition simCD, NormalModes simNM,
                                            CoordinateDefinition  otherCD, PotentialMasterList potentialMaster,
                                            int[] otherNCells, NormalModes otherNM) {
        this(sim, space, temp, simCD, simNM, otherCD, potentialMaster, 
                otherNCells, otherNM, "file");
    }
    
    public MeterDifferentImageAddDoubleSize(Simulation sim, Space space,
                                            double temp, CoordinateDefinition simCD, NormalModes simNM,
                                            CoordinateDefinition otherCD, PotentialMasterList potentialMaster,
                                            int[] otherNCells, NormalModes otherNM, String otherFilename){
        
        super(sim, space, temp, simCD, simNM, otherCD, potentialMaster, 
                otherNCells, otherNM, otherFilename);
    }
    
    public double getDataAsScalar() {
        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell[] cells = cDef.getBasisCells();
        double normalization = 1/Math.sqrt(2.0);  //DS!
        BasisCell cell = simCells[0];
        newU = new double[cDim];
        //Calculate the positions due to the simulation system.
        double[] simZero = simCDef.calcU(simCDef.getBasisCells()[0].molecules);
        
        //Create the last normal mode coordinates from the Gaussian distribution    
        for(int i = 0; i < gaussCoord.length; i++){
            gaussCoord[i] = random.nextGaussian() * sqrtTemperature;
        }
        
        //Calculate the positions due to the Gaussian normal mode coordinate
        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            int etaCount = 0;   //etaCount counts through "wv" for etas.
            //Calculate the change in positions.
            double kR = waveVectors[1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            if(wvCoeff[1] == 1.0){
                for (int iMode = 0; iMode < cDim; iMode++){
                    if(!(oneOverSqrtOmega2[1][iMode] == 0.0)){
                        for (int iCD = 0; iCD < cDim; iCD++){
                            newU[iCD] += sqrtWVC[1] * eigenVectors[1][iMode][iCD] 
                                * oneOverSqrtOmega2[1][iMode]
                                * (gaussCoord[etaCount] * coskR - gaussCoord[etaCount+1]* sinkR);
                        }
                        etaCount += 2;
                    }
                }
            } else {
                for (int iMode = 0; iMode < cDim; iMode++){
                    if(!(oneOverSqrtOmega2[1][iMode] == 0.0)){
                        for (int iCD = 0; iCD < cDim; iCD++){
                            newU[iCD] += sqrtWVC[1] * eigenVectors[1][iMode][iCD] 
                                * oneOverSqrtOmega2[1][iMode]
                                * (gaussCoord[etaCount] * coskR);
                        }
                        etaCount += 1;
                    }
                }
            }
            
            for (int i=0; i<cDim; i++) {
                newU[i] += simZero[i];
                newU[i] *= normalization;
            }
            
            cDef.setToU(cells[iCell].molecules, newU);
        }
        return meterPE.getDataAsScalar();
    }
    
}
