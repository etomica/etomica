package etomica;

import etomica.units.*;
import etomica.utility.Histogram;
import etomica.utility.History;
import etomica.utility.Function;

/**
 * A collection of Meters, or something that behaves as if it were a collection of Meters.
 * Exists to permit several properties to be measured at the same time, perhaps because
 * it is more efficient than computing them in separate meters.  Differs from a simple
 * Meter in that the data-access methods (currentValue, average, etc.), take an integer
 * argument indicating "which Meter's" value is requested.  Often there won't actually
 * exist a separate Meter for each measured property, but "PseudoMeters" can be constructed
 * and returned as an array of all meters by the allMeters() method.  Each meter in this
 * array behaves as a stand-alone meter would.
 *
 * @author David Kofke
 */
public abstract class MeterGroup extends MeterAbstract implements DataSource  {
    
    public static final String VERSION = "MeterGroup:01.05.24/"+MeterAbstract.VERSION;
    
    MeterAbstract.Accumulator[] accumulator;
    MeterScalar[] meters;
    protected double[] currentValues;
    private double[] values;
    private Function function;
    private int nMeters;
    protected String[] labels;
    
    private static final String[] sourcesAsText = new String[] {"History", "Histogram"};

	public MeterGroup(Simulation sim, int nMeters) {
	    super(sim);
	    setActive(true);  //default is to have meter do averages after some number of integrationIntervalEvents
        this.nMeters = nMeters;
        accumulator = new MeterAbstract.Accumulator[nMeters];
        for(int i=0; i<nMeters; i++) accumulator[i] = new MeterAbstract.Accumulator();
        values = new double[nMeters];
        currentValues = new double[nMeters];
        labels = new String[nMeters];
	}
	
	/**
	 * Method called to update the current values of the meters.
	 */
	 //a call to this method should leave the currentValues[] array filled with
	 //the latest values of the properties
	 public abstract void updateValues();
	 
	/**
	 * Number of meters in this group.
	 */
	public int nMeters() {return nMeters;}
	
	/**
	 * Sets a function that is applied to each meter's average.
	 * The meter's average() method returns the value returned by
	 * the function when applied to the value defined for measurement by the meter.
	 * For example, the function might square the meter's value, or divide by a temperature.
	 * If no function is specified (or is set to null), the meter returns its average, unaltered.
	 */
	public void setFunction(Function f) {function = f;}
	/**
	 * Returns the function that has been specified to apply to the meter's average.
	 */
	public Function getFunction() {return function;}
	
	/**
	 * Returns the current value of the i-th meter, as determined during
	 * the last call to updateValues.
	 * Note that the value returned by this method is not transformed by the function associated
	 * with this meter.
	 */
	public double currentValue(int i) {
	    return currentValues[i];
	}
	
    /**
     * Method to update running sums for averages and error statistics.
     * Called by integrationIntervalAction every updateInterval times an integration event is received
     */
	public void updateSums() {
	    updateValues();
	    for(int i=0; i<nMeters(); i++) {
	        accumulator[i].add(currentValues[i]);
	    }
	}
	
	/**
	 * For internal use.
	 */
	 protected Accumulator[] allAccumulators() {
	    return accumulator;
	 }
	
	public MeterScalar[] allMeters() {
	    if(meters == null) {
	        meters = new MeterScalar[nMeters];
	        for(int i=0; i<nMeters; i++) {
	            meters[i] = new PseudoMeter(simulation(), i);
	        }
	    }
	    return meters;
	}
	
	public double[] values(DataSource.ValueType type) {return values((MeterAbstract.ValueType)type);}
	/**
	 * Returns the value indicated by the argument.
	 */
	public double[] values(MeterAbstract.ValueType type) {
	    if(type==MeterAbstract.AVERAGE) 
	        for(int i=0; i<nMeters; i++) values[i] =  average(i);
	    else if(type==MeterAbstract.MOST_RECENT) 
	        for(int i=0; i<nMeters; i++) values[i] =  mostRecent(i);
	    else if(type==MeterAbstract.CURRENT) 
	        for(int i=0; i<nMeters; i++) values[i] =  currentValue(i);
	    else if(type==MeterAbstract.MOST_RECENT_BLOCK) 
	        for(int i=0; i<nMeters; i++) values[i] =  mostRecentBlock(i);
	    else if(type==MeterAbstract.ERROR) 
	        for(int i=0; i<nMeters; i++) values[i] =  error(i);
	    else if(type==MeterAbstract.VARIANCE) 
	        for(int i=0; i<nMeters; i++) values[i] =  variance(i);
	    return values;
	}
	
	public String[] getSourcesAsText() {
	    int nSources = 0;
	    if(historying) nSources++;
	    if(histogramming) nSources++;
	    String[] sources = new String[nSources];
	    int i=0;
	    if(historying) sources[i++] = "History";
	    if(histogramming) sources[i++] = "Histogram";
	    return sources;
	}
	public DataSource[] getDataSources(String text) {
	    DataSource[] sources = new DataSource[nMeters];
	    if(text.equals("History")) 
	        for(int i=0; i<nMeters; i++) {
	            sources[i] = accumulator[i].history();
	            ((History)sources[i]).setLabel(labels[i]);
	        }
	    else if(text.equals("Histogram")) 
	        for(int i=0; i<nMeters; i++) {
	            sources[i] = accumulator[i].histogram();
	            ((Histogram)sources[i]).setLabel(labels[i]);
	        }
	    return sources;
	}
	
	public History getHistory(int i) {
	    return accumulator[i].history();
	}
	
    /**
     * Returns the current value of the average, transformed by this meter's function, if defined.
     */
	public double average(int i) {
	    return (function==null) ? accumulator[i].average() : function.f(accumulator[i].average());
	}

    /**
     * Returns the current value of the variance.
     * No function transformation is applied.
     */
	public double variance(int i) {return accumulator[i].variance();}
	
	/**
	 * Returns the current value of the error bar (67% confidence limit),
	 * appropriately transformed by the meter's function, if defined.
	 */
	public double error(int i) {
	    if(function == null) return accumulator[i].error();
	    else {//have not carefully considered if this is correct
	        return Math.abs(function.dfdx(accumulator[i].average()))*accumulator[i].error();
	    }
	}
	
	/**
	 * Returns the value passed most recently obtained via an integrator intervalEvent
	 * (does not give the value last returned by any direct call to currentValue
	 *  i.e. by a call to currentValue that was not made in response to an interval event).
	 * Transformed by meter's function, if defined.
	 */
	public double mostRecent(int i) {
	    return (function==null) ? accumulator[i].mostRecent() : function.f(accumulator[i].mostRecent());}
	
    /**
     * Returns the most recent block average, or NaN if no block averages are complete.
     * Transformed by the meter's function, if defined.
     */
     public double mostRecentBlock(int i) {
        return (function==null) ? accumulator[i].mostRecentBlock() : function.f(accumulator[i].mostRecentBlock());}
	
	
	/**
	 * Returns the current histogram of measured values.
	 * Histogram is recorded only if histogramming is set to true.  
	 * If histogram was never kept for this meter, an all-zero histogram is returned.
	 */
	 public Histogram histogram(int i) {return accumulator[i].histogram();}
	 
	 /**
	  * Meter facade that gives the impression of providing data independently,
	  * although it is actually serving as a wrapper for the one of the data values
	  * collected by the meter group.
	  */
	 public class PseudoMeter extends MeterScalar {
       private final int index;
        PseudoMeter(Simulation sim, int i) {
            super(sim);
            index = i;
            this.accumulator = MeterGroup.this.accumulator[i];
        } 
        public double  currentValue() {return MeterGroup.this.currentValue(index);}
        public boolean usesPhaseIteratorFactory() {return false;}
        public boolean usesPhaseBoundary() {return false;}
        public String getLabel() {return labels[index];}
        public Dimension getDimension() {return Dimension.NULL;}//temporary
        public void updateSums() {/*do nothing since this is taken care of by the group*/}
	 }//end of PseudoMeter
	 
}//end of MeterGroup class	 
