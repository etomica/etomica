package etomica.models.nitrogen;

import java.io.Serializable;

import etomica.action.IAction;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.HistogramExpanding;

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

    public MeterNormalizedCoord(IBox newBox, CoordinateDefinition coordDef, ISpecies species) {
        tag = new DataTag();
        this.coordinateDefinition = coordDef;
        
        if (species==null){
        	throw new RuntimeException("Need to set species in MeterNormalizedCoord class");
        }
        dof = coordDef.getCoordinateDim() /newBox.getNMolecules(species);
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
                                  
	        for (int i=0; i<dof; i++){ 
	        	for (int j=0; j<numMolecules; j++){
	            	histogramU[i].addValue(u[i+dof*j]);
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

    
}
