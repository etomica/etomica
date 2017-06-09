/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Space;
import etomica.units.Length;
import etomica.units.Null;
import etomica.data.histogram.HistogramExpanding;

/**
 * Calculates the average atomic displacement from their lattice sites
 */
public class MeterAtomicDisplacement implements IEtomicaDataSource, DataSourceIndependent, IAction {

    public MeterAtomicDisplacement(Space space, CoordinateDefinition coordinateDefinition) {
    	
    	this.coordinateDefinition = coordinateDefinition;
    	this.space = space;
    
    	histogram = new HistogramExpanding(0.002);
    	dataInfoFunction = new DataInfoFunction("Probability", Null.DIMENSION, this);
        tag = new DataTag();

        
    	distanceR = new DataDoubleArray(0);
        rDataInfo = new DataInfoDoubleArray("Atomic Distance", Length.DIMENSION, new int[]{0});
        rTag = new DataTag();
        rDataInfo.addTag(rTag);
    }
    	


    public Box getBox() {
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
    			IAtom a = molecules.getMolecule(i).getChildList().getAtom(0);
    			Vector pos = a.getPosition();
    			Vector site = coordinateDefinition.getLatticePosition(a);
    			
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
		
        //if (rDataInfo != null || rDataInfo.getLength() != histogram.getNBins()) {
            distanceR = new DataDoubleArray(histogram.getNBins());
           // rDataInfo = new DataInfoDoubleArray("Atomic Distance", Length.DIMENSION, new int[]{histogram.getNBins()});
          //  rDataInfo.addTag(tag);
        //}
		distanceR.E(histogram.xValues());
		
		return distanceR;
	}



	public DataInfoDoubleArray getIndependentDataInfo(int i) {
	    //if (rDataInfo != null || rDataInfo.getLength() != histogram.getNBins()) {
           // distanceR = new DataDoubleArray(histogram.getNBins());
	        rDataInfo = new DataInfoDoubleArray("Atomic Distance", Length.DIMENSION, new int[]{histogram.getNBins()});
	      //  rDataInfo.addTag(tag);
	    //}
		return rDataInfo;
	}

	public DataTag getIndependentTag() {
	    return rTag;
	}
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final DataTag tag, rTag;
    private DataInfoFunction dataInfoFunction;
    private DataInfoDoubleArray rDataInfo;
    protected Vector workVector;
    private DataFunction data;
    private Space space;
    private HistogramExpanding histogram;
    private DataDoubleArray distanceR;
	

}
