package etomica.virial.overlap;

import etomica.MeterScalar;
import etomica.SimulationElement;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * Meter to computer cluster integral by processing averages from overlap-
 * sampling meters.
 */
public class MeterClusterOverlap extends MeterScalar {

	/**
	 * Constructor for MeterClusterOverlap.
	 * @param parent
	 */
	public MeterClusterOverlap(SimulationElement parent, 
			MeterOverlapReference referenceMeter,
			MeterOverlapTarget targetMeterPositive, 
			MeterOverlapTarget targetMeterNegative) {
		super(parent);
		this.referenceMeter = referenceMeter;
		this.targetMeterPositive = targetMeterPositive;
		this.targetMeterNegative = targetMeterNegative;
	}

	/**
	 * @see etomica.MeterScalar#currentValue()
	 */
	public double currentValue() {
		double[] refNegative = referenceMeter.allMeters()[0].average();
		double[] refPositive = referenceMeter.allMeters()[1].average();
		double[] targetNegative = targetMeterNegative.average();
		double[] targetPositive = targetMeterPositive.average();
		double[] x = targetMeterPositive.X();
		double xN = intersection(x, refNegative, targetNegative);
		double xP = intersection(x, refPositive, targetPositive);
		double value =  Math.exp(xP) - Math.exp(xN);
//		System.out.println(value);
		return value;
	}
	
	private double intersection(double[] x, double[] y1, double[] y2) {
		//look for intersection as point where y1 goes from above(below) y2 to
		//below(above) y2
		boolean greater = y1[0] > y2[0];
		for(int i=1; i<x.length; i++) {
			if( (y1[i]>y2[i]) != greater) {//found it
				//linear interpolation for point of intersection
				double xI = x[i-1] + (y2[i]-y1[i])*(x[i]-x[i-1])/((y1[i]-y1[i-1])-(y2[i]-y2[i-1]));
				return xI;
			}
		}
		//did not find intersection
		return Double.NaN;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return null;
	}

	MeterOverlapReference referenceMeter;
	MeterOverlapTarget targetMeterPositive, targetMeterNegative;
}
