/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;
import etomica.virial.cluster.ClusterAbstract;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirial implements IDataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters) {
		clusters = aClusters;
        data = new DataDoubleArray(clusters.length);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{clusters.length});
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
        int j =0;
      //  System.out.println(pi);
        if (pi == 0 || pi == Double.POSITIVE_INFINITY || Double.isNaN(pi)) throw new RuntimeException("oops "+pi+" for box "+box.getIndex());
        double x[] = data.getData();
        for (int i=0; i<clusters.length; i++) {
            x[i] = clusters[i].value(box)/pi;
           // System.out.println(x[i] + " " + j);
            j++;
            if (Double.isNaN(x[i]) || Double.isInfinite(x[i])) throw new RuntimeException("oops for box "+box.getIndex()+" "+clusters[i].value(box)+" "+x[i]+" "+pi);
        }
//        System.out.println(box.getIndex()+" "+pi+" "+x[0]+" "+x[1]);
        return data;
    }
    
    public ClusterAbstract[] getClusters() {
        return clusters;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final ClusterAbstract clusters[];
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
