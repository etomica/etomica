/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.IAction;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Null;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

/**
 * Meter that measures the spherical coordinate angles
 * angle[0]: is theta, angle between the projection vector (molecular axis vector onto xy plane)
 * 				and x-axis
 * angle[1]: is phi,   angle between the molecular axis vector and z-axis  
 * 
 * @author Tai Boon Tan
 *
 */
public class MeterRotationDistributionGroup implements IDataSource, IAction, Serializable {

    public MeterRotationDistributionGroup(Box newBox, CoordinateDefinitionNitrogen coordDef) {
        tag = new DataTag();
        this.coordinateDefinition = coordDef;
 
        molAxis = Space3D.makeVector(3);
        initOrient = Space3D.makeVector(3);
        
        uDistributions = new DataDoubleArray[2];
        uDistributionsInfo = new DataInfoDoubleArray[2];
        
        histogramNotSimple = new HistogramNotSoSimple(new DoubleRange (0, Math.PI/2.0));
        histogramNotSimple.setDoAveraging(false);
        histogramExpanding = new HistogramExpanding(0.0005);
        
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Assigning the histogram to data
     */
    public void actionPerformed() {
    
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        
        double cosThetaTotal = 0.0;
        for (int iCell = 0; iCell<cells.length; iCell++) {
      
            BasisCell cell = cells[iCell];
            IMoleculeList molecules = cell.molecules;
            int numMolecules = molecules.size();
            
            for (int iMol=0; iMol<numMolecules; iMol++){
            	
	          	IMolecule molecule = molecules.get(iMol);
	          	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
		    	Vector leafPos1 = molecule.getChildList().get(1).getPosition();
		
		    	molAxis.Ev1Mv2(leafPos1, leafPos0);
		       	molAxis.normalize();
		       	
		       	initOrient.E(coordinateDefinition.getMoleculeOrientation(molecule)[0]);
		       	double cosT = molAxis.dot(initOrient);
	            angle = Math.acos(cosT);
		       	
//	            System.out.println(angle + " "+Math.sin(angle) + " " +1/Math.sin(angle) + " " + angle);
	            histogramNotSimple.addValue(angle, 1/Math.sin(angle));
	            
	            cosThetaTotal += cosT;
            }
        } //end of cell; there is only 1 cell
        
        
        histogramExpanding.addValue(cosThetaTotal/(cells.length*cells[0].molecules.size()));
    }

    /**
     * Returns the DataGroup of uDistribtion[i]
     */
    public IData getData() {
 
		int nBinNotSimple = histogramNotSimple.getNBins();
		uDistributions[0] = new DataDoubleArray(new int[]{2, nBinNotSimple});
        uDistributionsInfo[0] = new DataInfoDoubleArray("Angle", Null.DIMENSION, new int[]{2, nBinNotSimple});
        
		uDistributions[0].assignColumnFrom(0, histogramNotSimple.xValues());
    	uDistributions[0].assignColumnFrom(1, histogramNotSimple.getHistogram());
    		    	
    	
    	
    	int nBinExpanding = histogramExpanding.getNBins();
    	uDistributions[1] = new DataDoubleArray(new int[]{2, nBinExpanding});
        uDistributionsInfo[1] = new DataInfoDoubleArray("Angle", Null.DIMENSION, new int[]{2, nBinExpanding});
        
		uDistributions[1].assignColumnFrom(0, histogramExpanding.xValues());
    	uDistributions[1].assignColumnFrom(1, histogramExpanding.getHistogram());
    	
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
		String fNameRot = filename+".RotDistrb";
		String fNameCons = filename+".ConsDistrb";
		try {
			FileWriter fileWriterRot = new FileWriter(fNameRot,false);
			FileWriter fileWriterCons = new FileWriter(fNameCons,false);
			
			DataDoubleArray uDistribution0 = (DataDoubleArray)uData.getData(0);
			DataDoubleArray uDistribution1 = (DataDoubleArray)uData.getData(1);
			
			
			for (int j=0; j<uDistribution0.getLength()/uDistribution0.getArrayDimension(); j++){
				fileWriterRot.write(uDistribution0.getValue(new int[]{0,j})+" "+ uDistribution0.getValue(new int[]{1,j}) + "\n");
			}
		
			for (int j=0; j<uDistribution1.getLength()/uDistribution1.getArrayDimension(); j++){
				fileWriterCons.write(uDistribution1.getValue(new int[]{0,j})+" "+ uDistribution1.getValue(new int[]{1,j}) + "\n");
			}
			
			fileWriterRot.close();
			fileWriterCons.close();
			
		} catch(IOException e){
			throw new RuntimeException("<MeterRotationDistributionGroup> Failed to write ratation distribution" + e);
		
		}
		
	}
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinitionNitrogen coordinateDefinition;
    private final DataTag tag;
    private IDataInfo dataInfo;
    private DataGroup data;
    private HistogramNotSoSimple histogramNotSimple;
    private HistogramExpanding histogramExpanding;
    
    
    private DataDoubleArray[] uDistributions;
    private DataInfoDoubleArray[] uDistributionsInfo;
    private Vector initOrient;
    private Vector molAxis;
    private double angle;
    
}
