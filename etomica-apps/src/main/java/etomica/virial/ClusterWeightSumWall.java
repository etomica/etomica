/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.data.IData;
import etomica.data.IDataSource;

public class ClusterWeightSumWall implements ClusterWeight {
	
	public ClusterWeightSumWall(DataSourceClusterWall meter, int pointCount ) {
		this.meter = meter; 
		this.pointCount = pointCount;
		
	}
	
	public ClusterAbstract makeCopy() {
		
		return new ClusterWeightSumWall(meter, pointCount);
	}

	public int pointCount() {
		
		return pointCount;
	}

	public double value(BoxCluster box) {
		meter.setBox(box);
		IData idata = meter.getData();
		return idata.getValue(idata.getLength()-1);
	}

	public void setTemperature(double temperature) {
		
		
	}
	private final DataSourceClusterWall meter;
	private final int pointCount;

    public interface DataSourceClusterWall extends IDataSource {
        public void setBox(BoxCluster box);
	}
	
}
