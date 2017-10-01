/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterSamplingWeight implements IDataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterSamplingWeight(DataSourceAtomDistance dataDistance) {
		this.dataDistance = dataDistance; 
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
        double pi = box.getSampleCluster().value(box);
        data.getData()[0] = dataDistance.getDataAsScalar();
        data.getData()[1] = pi;
        return data;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final DataSourceAtomDistance dataDistance;
	private final DataDoubleArray data;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
