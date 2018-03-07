/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialExternalFieldConfined implements IDataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialExternalFieldConfined(ClusterAbstract aCluster, double[] wallposition, double [] walldistance) {
		cluster = aCluster;
		this.wallPosition = wallposition;
		wallDistance = walldistance;
        data = new DataDoubleArray(new int[]{walldistance.length+1,wallposition.length+1});
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{walldistance.length+1,wallposition.length+1});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}

	public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
        double pi = box.getSampleCluster().value(box);
        double x[] = data.getData();
        IAtomList atoms = box.getLeafList();
        double lowestatom = 0;
        double highestatom = 0;
        for (int i = 1; i<atoms.size(); i++){
        	double z = atoms.get(i).getPosition().getX(2);
        	if (z<lowestatom){
        		lowestatom = z;      		
        	}
        	else if (z>highestatom){
        		highestatom = z;       		
        	}
        }
        
        double v=cluster.value(box)/pi;
        int n=wallDistance.length+1;
        for (int i=0; i<wallPosition.length; i++) {
        	x[i*n]=v;
            if (lowestatom-0.5 < wallPosition[i]){
            	x[i*n] -= 0.5*v;
            }
            if (-highestatom-0.5 < wallPosition[i]){
            	x[i*n] -= 0.5*v; 
            }
            for (int j=0; j<wallDistance.length; j++){
            	x[i*n+j+1]=v;
            	if (lowestatom-0.5 < wallPosition[i]||highestatom+0.5 > wallPosition[i]+wallDistance[j]){
            		x[i*n+j+1] = 0;
            	}
            }
            
        }
        
        for (int j=0; j<wallDistance.length; j++){
        	double l=(wallDistance[j]-1)-(-lowestatom+highestatom);
        	x[n*wallPosition.length+j+1]=l<0 ? 0 : l*v;

        }
        x[n*wallPosition.length]=(atoms.size()-1+0.5*(lowestatom-highestatom))*v;
        
        return data;
    }
    
    public ClusterAbstract getCluster() {
        return cluster;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final ClusterAbstract cluster;
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private final double [] wallPosition;
    private final double [] wallDistance;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
