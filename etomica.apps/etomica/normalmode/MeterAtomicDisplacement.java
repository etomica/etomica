package etomica.normalmode;

import etomica.action.IAction;
import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.ISpace;
import etomica.units.Length;
import etomica.units.Null;
import etomica.util.Histogram;
import etomica.util.HistogramExpanding;

/**
 * Calculates the average atomic displacement from their lattice sites
 */
public class MeterAtomicDisplacement implements IEtomicaDataSource, DataSourceIndependent, IAction {

    public MeterAtomicDisplacement(ISpace space, CoordinateDefinition coordinateDefinition) {
    	
    	this.coordinateDefinition = coordinateDefinition;
    	this.space = space;
    
    	histogram = new HistogramExpanding(0.002);
    	dataInfoFunction = new DataInfoFunction("Probability", Null.DIMENSION, this);
        tag = new DataTag();

        
    	distanceR = new DataDoubleArray(getBox().getLeafList().getAtomCount());
    }
    	


    public IBox getBox() {
        return coordinateDefinition.getBox();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfoFunction;
    }

    /**
     * 
     */
    public void actionPerformed() {
    	
    	for (int iCell = 0; iCell < coordinateDefinition.cells.length; iCell++){
    		BasisCell cell = coordinateDefinition.getBasisCells()[iCell];
    		IMoleculeList molecules = cell.molecules;
    		
    		
    		for (int i=0; i< cell.molecules.getMoleculeCount(); i++){
    			workVector = space.makeVector();
    			IAtomPositioned a = (IAtomPositioned)molecules.getMolecule(i).getChildList().getAtom(0);
    			IVectorMutable pos = a.getPosition();
    			IVectorMutable site = coordinateDefinition.getLatticePosition((IAtom)a);
    			
    			workVector.Ev1Mv2(pos, site);
    			
    			histogram.addValue(Math.sqrt(workVector.squared()));
    		}
    	}
    	
    }
    
    public IData getData(){
    	double[] rData = new double[histogram.getNBins()];
    	data = new DataFunction(new int[]{histogram.getNBins()});
    	//dx is the half the bin size
    	/*
    	 * similar to the implementation in MeterRDF
    	 * normalizing the probability
    	 */
    	double dx = histogram.getDeltaX()/2;
    	double sum = 0;
    	for (int i=0; i<histogram.getNBins(); i++){
    		rData[i] = histogram.getHistogram()[i];
    		double vShell = space.sphereVolume(histogram.xValues()[i]+dx) -space.sphereVolume(histogram.xValues()[i]-dx);
    		rData[i] /= vShell;
    		sum += rData[i]*dx*2;
    	}
    	for (int i=0; i<histogram.getNBins(); i++){
    		rData[i] /= sum;
    	}
    	
    	data.E(rData);
    	
        return data;
    }

   
    public DataTag getTag() {
        return tag;
    }
    
    
	public int getIndependentArrayDimension() {
		return 1;
	}



	public DataDoubleArray getIndependentData(int i) {
		
		distanceR = new DataDoubleArray(histogram.getNBins());
		distanceR.E(histogram.xValues());
		
		return distanceR;
	}



	public DataInfoDoubleArray getIndependentDataInfo(int i) {
		
		dataInfo = new DataInfoDoubleArray("Atomic Distance", Length.DIMENSION, new int[]{histogram.getNBins()});
		return dataInfo;
	}
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final DataTag tag;
    private DataInfoFunction dataInfoFunction;
    private DataInfoDoubleArray dataInfo;
    protected IVectorMutable workVector;
    private DataFunction data;
    private ISpace space;
    private HistogramExpanding histogram;
    private DataDoubleArray distanceR;
	

}
