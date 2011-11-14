package etomica.virial;

import etomica.data.IData;

public class ClusterWeightSumWall implements ClusterWeight {
	
	public ClusterWeightSumWall(MeterVirialExternalFieldOverlapRho meter, int pointCount ) {
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
		double sumwall = 0;
		for (int i=1; i<idata.getLength()-1; i++) {
			sumwall +=Math.abs(idata.getValue(i));
		}
		
		return sumwall;
	}

	public void setTemperature(double temperature) {
		
		
	}
	private final MeterVirialExternalFieldOverlapRho meter;
	private final int pointCount;
}
