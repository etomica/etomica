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
public class MeterVirialExternalField implements IDataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialExternalField(ClusterAbstract aCluster, double[] wallposition) {
		cluster = aCluster;
		this.wallPosition = wallposition;
        data = new DataDoubleArray(wallposition.length+1);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{wallposition.length+1});
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
        for (int i = 1; i<atoms.size(); i++){
        	double z = atoms.get(i).getPosition().getX(2);
        	if (z<lowestatom){
        		lowestatom = z;
        		
        	}
        }
        double v=cluster.value(box)/pi;
        for (int i=0; i<wallPosition.length; i++) {
            if (lowestatom-0.5 < wallPosition[i]){
            	x[i] = 0;
            }
            else x[i]=v;
        }
        x[x.length-1]=(atoms.size()-1+lowestatom)*v;
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
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
