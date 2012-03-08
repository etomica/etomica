package etomica.virial;

import etomica.data.IData;

public class ClusterWeightSumWall implements ClusterWeight {
	
	public ClusterWeightSumWall(MeterVirialExternalFieldOverlapConfined meter, int pointCount ) {
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
	private final MeterVirialExternalFieldOverlapConfined meter;
	private final int pointCount;
}
