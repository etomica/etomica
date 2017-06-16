/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.IAction;
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
import etomica.species.ISpecies;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;

import java.io.Serializable;

/**
 * Sample the distribution of the normalized coordinate for the atoms/ molecules
 *  from coordinate definition, calcU method. They are the quantities of deviation
 *  in translation and rotation from the atomic/molecular nominal position and
 *  orientation
 *  
 *  Eg. N2 molecule with 5 D.O.F, this class will return 5 distributions of 
 *  	3 translation and 2 rotation motion respectively.
 */
public class MeterNormalizedCoord implements IEtomicaDataSource, IAction, Serializable {
	 public MeterNormalizedCoord(Box newBox, CoordinateDefinition coordDef, ISpecies species) {
		 this(newBox, coordDef, species,false);
	 }
	
    public MeterNormalizedCoord(Box newBox, CoordinateDefinition coordDef, ISpecies species, boolean isVolFluctuation) {
        tag = new DataTag();
        this.coordinateDefinition = coordDef;
        this.isVolFluctuation = isVolFluctuation;
        
        if (species==null){
        	throw new RuntimeException("Need to set species in MeterNormalizedCoord class");
        }
        /*
         * minus 3 is to ignore the volume fluctuation*****!!
         */
        if(isVolFluctuation){
        	dof = (coordDef.getCoordinateDim()- 1) /newBox.getNMolecules(species) +1;
        } else {
        	dof = coordDef.getCoordinateDim() /newBox.getNMolecules(species);
            
        }
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
            
            double[] u = coordinateDefinition.calcU(molecules);
            
            if(isVolFluctuation){
            	for (int i=0; i<(dof-1); i++){ 
 		        	for (int j=0; j<numMolecules; j++){
 		            	histogramU[i].addValue(u[i+(dof-1)*j]);
 		            }
 		        } 
            	
            	//for (int i=0; i<3; i++){
            		histogramU[(dof-1)].addValue(u[(coordinateDefinition.getCoordinateDim()-1)]);
            	//}
            	
            } else {
		        for (int i=0; i<dof; i++){ 
		        	for (int j=0; j<numMolecules; j++){
		            	histogramU[i].addValue(u[i+dof*j]);
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
    
	public boolean isVolFluctuation() {
		return isVolFluctuation;
	}

	public void setVolFluctuation(boolean isVolFluctuation) {
		this.isVolFluctuation = isVolFluctuation;
	}
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private String name;
    private final DataTag tag;
    private IEtomicaDataInfo dataInfo;
    private DataGroup data;
    private HistogramExpanding[] histogramU;
    private DataDoubleArray[] uDistributions;
    private DataInfoDoubleArray[] uDistributionsInfo;
    private int dof;
    private boolean isVolFluctuation;


    
}
