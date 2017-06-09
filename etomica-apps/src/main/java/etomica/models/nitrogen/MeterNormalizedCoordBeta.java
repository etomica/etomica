/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.IAction;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

/**
 * Sample the distribution of the normalized coordinate for the atoms/ molecules
 *  from coordinate definition, calcU method. They are the quantities of deviation
 *  in translation and rotation from the atomic/molecular nominal position and
 *  orientation
 *  
 *  Specific class to sample the distribution of the translational dof of beta-N2
 *   with 2 molecules in a unit cell and there are 2 distinct unit cells. So the
 *   total number of parameter that we sample for the system is 3 x 2 x 2 = 12.
 *  
 */
public class MeterNormalizedCoordBeta implements IEtomicaDataSource, IAction, Serializable {

    public MeterNormalizedCoordBeta(Box newBox, CoordinateDefinition coordDef, ISpecies species) {
        tag = new DataTag();
        this.coordinateDefinition = coordDef;
        
        if (species==null){
        	throw new RuntimeException("Need to set species in MeterNormalizedCoord class");
        }
        dof = 12;
         
        uDistributions = new DataDoubleArray[dof];
        uDistributionsInfo = new DataInfoDoubleArray[dof];
        
        histogramU = new HistogramExpanding[dof];
        for (int i=0; i<dof; i++){
        	histogramU[i] = new HistogramExpanding(0.005);
        }
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Assigning the histogram to data
     */
    public void actionPerformed() {

        BasisCell[] cells = coordinateDefinition.getBasisCells();
                
        for (int iCell = 0; iCell<cells.length; iCell++) {
        	
            BasisCell cell = cells[iCell];
            IMoleculeList molecules = cell.molecules;
            int numMolecules = molecules.getMoleculeCount();
            int nCell = (int)Math.pow(numMolecules/1.999999999, 1.0/3.0);
            double[] u = coordinateDefinition.calcU(molecules);
            
            //unit cell A
	        for (int i=0; i<6; i++){ 
			    if(i<3){	
			    	for (int j=0; j<numMolecules; j+=2){
			    		if(j>0 && j%(nCell*2)==0){
			    			j+=(nCell*2);
			    			if(j>=numMolecules){
								break;
							}
			    		}
			    		histogramU[i].addValue(u[j*5+i]);
			        }
			    }else{
			    	for (int j=1; j<numMolecules; j+=2){
			    		if(j>1 && j%(nCell*2)==1){
			    			j+=(nCell*2);
			    			if(j>=numMolecules){
								break;
							}
			    		}
			    		histogramU[i].addValue(u[j*5+(i-3)]);
			        }
			    }
		    } 
			
	        //unit cell B
		    for (int i=6; i<12; i++){ 
			    if(i<9){	
			    	for (int j=(nCell*2); j<numMolecules; j+=2){
			    		if(j>(nCell*2) && j%(nCell*2)==0){
			    			j+=(nCell*2);
			    			if(j>=numMolecules){
								break;
							}
			    		}
			    		histogramU[i].addValue(u[j*5+(i-6)]);
			        }
			    }else{
			    	for (int j=((nCell*2)+1); j<numMolecules; j+=2){
			    		if(j>((nCell*2)+1) && j%(nCell*2)==1){
			    			j+=(nCell*2);
			    			if(j>=numMolecules){
								break;
							}
			    		}
			    		histogramU[i].addValue(u[j*5+(i-9)]);
			        }
			    }
		    } 
	        
        } //end of cell; there is only 1 cell
    }

    /**
     * Returns the DataGroup of uDistribtion[i]
     */
    public IData getData() {
 
    	CompoundDimension length = new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{1});
        
    	for(int i=0; i<dof; i++){
    		int nBin = histogramU[i].getNBins();
    		uDistributions[i] = new DataDoubleArray(new int[] {2, nBin});
    		uDistributionsInfo[i] = new DataInfoDoubleArray("u", length, new int[]{2, nBin});
    		
    		uDistributions[i].assignColumnFrom(0, histogramU[i].xValues());
    		uDistributions[i].assignColumnFrom(1, histogramU[i].getHistogram());
    		
    	}
    	
    	data = new DataGroup(uDistributions);
        dataInfo = new DataInfoGroup("all uDistributions", Null.DIMENSION, uDistributionsInfo);
      
    	return data;
    	
    }

    /**
     * Sets the tensor summation to 0.
     */
    public void reset() {
        data.E(0);
    }

    public DataTag getTag() {
        return tag;
    }
    
	public void writeUdistribution(String filename){
		DataGroup uData = (DataGroup)getData();
		
		for (int i=0; i<uData.getNData(); i++){
			String fName = filename+"U"+i+".coord";
			try {
				FileWriter fileWriter = new FileWriter(fName,false);
				
				DataDoubleArray uDistribution = (DataDoubleArray)uData.getData(i);
				
				for (int j=0; j<uDistribution.getLength()/uDistribution.getArrayDimension(); j++){
					fileWriter.write(uDistribution.getValue(new int[]{0,j})+" "+ uDistribution.getValue(new int[]{1,j}) + "\n");
				}
			
				fileWriter.close();
				
			} catch(IOException e){
				throw new RuntimeException("Failed to write coord data orientation U" + e);
			
			}
		}
		
	}
    
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final DataTag tag;
    private IEtomicaDataInfo dataInfo;
    private DataGroup data;
    private HistogramExpanding[] histogramU;
    private DataDoubleArray[] uDistributions;
    private DataInfoDoubleArray[] uDistributionsInfo;
    private int dof;

    
}
