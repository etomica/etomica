/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Returns cluster and distribution function values from ClusterSumIC 
 */
public class MeterVirialIC implements IDataSource, java.io.Serializable {
    protected static int count;
    protected final ClusterSumIC targetCluster;
    private final DataDoubleArray data;
    private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialIC(ClusterSumIC targetCluster) {
	    this.targetCluster = targetCluster;
	    int nv = targetCluster.getNumValues();
        data = new DataDoubleArray(nv);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{nv});
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
        if (pi == 0 || pi == Double.POSITIVE_INFINITY || Double.isNaN(pi)) {
            box.trialNotify();
            box.getSampleCluster().value(box);
            throw new RuntimeException("oops "+pi);
        }
        double x[] = data.getData();
        int nv = targetCluster.getNumValues();
        x[nv-1] = targetCluster.value(box)/pi;
        for (int i=0; i<nv-1; i++) {
            x[i] = targetCluster.value12(i)/pi;
        }
        return data;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

}
