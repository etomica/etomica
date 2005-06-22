package etomica.data.meter;

import etomica.Accumulator;
import etomica.Meter;
import etomica.Simulation;
import etomica.data.DataSourceFunction;
import etomica.units.Dimension;

/**
 * A collection of MeterFunctions, or something that behaves as if it were a
 * collection of MeterFunctions. Exists to permit several properties to be
 * measured at the same time, perhaps because it is more efficient than
 * computing them in separate meters.  MeterFunction instances are accessess
 * through the allMeters method, which returns an array of MeterFunction
 * instances.
 *
 * @author David Kofke
 */
public abstract class MeterFunctionGroup extends Meter  {
    
	public static final String VERSION = "MeterFunctionGroup:03.07.21/"+MeterAbsMeter;
    
	protected PseudoMeter[] meters;
	protected int nMeters;
	protected double[][] y;
	protected int nPoints;
    
	public MeterFunctionGroup(int nMeters) {
		super();
		setActive(true);  //default is to have meter do averages after some number of integrationIntervalEvents
		this.nMeters = nMeters;
		makeMeters();
	}
	
	/**
	 * Method called to update the current values of the meters.
	 * A call to this method should leave the y[][] array filled with the latest
	 * values of the properties
	 */
	 public abstract void updateValues();
	 
	/**
	 * Number of meters in this group.
	 */
	public int nMeters() {return nMeters;}
		
	/**
	 * Method to update running sums for averages and error statistics.
	 * Called by integrationIntervalAction every updateInterval times an integration event is received
	 */
	public void updateSums() {
		updateValues();
		for(int i=0; i<nMeters(); i++) {
			meters[i].updateSums();
		}
	}
	
	/**
	 * For internal use.
	 */
	 protected Accumulator[] allAccumulators() {
	 	int count = 0;
	 	for(int i=0; i<nMeters; i++) count += meters[i].getNPoints();
	 	Accumulator[] accumulator = new Accumulator[count];
	 	int k = 0;
	 	for(int i=0; i<nMeters; i++) {
	 		int n = meters[i].getNPoints();
	 		Accumulator[] a = meters[i].accumulator();
	 		for(int j=0; j<n; j++) accumulator[k++] = a[j];
	 	}
		return accumulator;
	 }
	
	public DataSourceFunction[] allMeters() {
		return meters;
	}
	
	private void makeMeters() {
		meters = new PseudoMeter[nMeters];
		for(int i=0; i<nMeters; i++) {
			meters[i] = new PseudoMeter(simulation(), i);
		}
		updateYArrays();
	}
	
	private void updateYArrays() {
		y = new double[nMeters][];
		for(int i=0; i<nMeters; i++) {
			y[i] = ((PseudoMeter)meters[i]).y();
		}
	}
	
	/**
	 * Calls setX for all meters in group, using given arguments.
	 */
	public void setX(double min, double max, int n) {
		nPoints = n;
		for(int i=0; i<nMeters; i++) {meters[i].setX(min, max, n);}
	}

		
	 /**
	  * Meter facade that gives the impression of providing data independently,
	  * although it is actually serving as a wrapper for the one of the data values
	  * collected by the meter group.
	  */
	 public class PseudoMeter extends DataSourceFunction {
	   private final int index;
		PseudoMeter(Simulation sim, int i) {
			super(sim);
			index = i;
			setActive(false);
		}
		Accumulator[] accumulator() {return accumulator;}
		public double[] getData() {return this.y;}//should be preceded by Group's call to updateSums
		public Dimension getDimension() {return Dimension.NULL;}//temporary
		public Dimension getXDimension() {return Dimension.NULL;}//temporary
		protected void resizeArrays() {
			super.resizeArrays();
			if(MeterFunctionGroup.this != null) updateYArrays();//need to check for null because PseudoMeter is constructed in constructor of MeterFunctionGroup, but MeterFunction parent class invokes setX, which leads to this resizeArrays method being invoked before enclosing class is constructed
		}
		private double[] y() {return this.y;}
	 }//end of PseudoMeter
	 
}//end of MeterGroup class	 
