package etomica;

import etomica.utility.Histogram;
import etomica.utility.History;
import etomica.utility.Function;

/**
 * Meter for recording and averaging a simple scalar of type double.
 * Subclasses of this abstract class need only to define the currentValue method
 * to calculate the desired property.  All averaging and error analysis is handled automatically by 
 * the methods defined by the MeterAbstract class.
 * The default condition of a Meter is "active", meaning that it increments sums for averages upon
 * receiving sufficient interval events from the integrator.
 */
 
 /* History
  * 10/17/02 (DAK) added Molar inner class that defines a function that makes Meter
  *          return averages of molar property (instead of extensive value)
  */
public abstract class MeterScalar extends MeterAbstract implements DataSource.Wrapper, DatumSource  {
    
    public static final String VERSION = "Meter:01.05.14/"+MeterAbstract.VERSION;
    
    MeterAbstract.Accumulator accumulator = new MeterAbstract.Accumulator();
    private Function function;
    
	public MeterScalar(SimulationElement parent) {
	    super(parent);
	    setActive(true);  //default is to have meter do averages after some number of integrationIntervalEvents
	}
	
	/**
	 * Sets a function that is applied to the meter's average.
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
	 * Defined by the subclass to specify what property is measured by the meter.
	 * The method is called by the updateSums method to add the current value to the accumulator.
	 * This method may be called directly by any object to obtain the current value of the measured
	 * property.  As an example, a energy meter is present in each phase for the purpose of providing
	 * a value for the potential energy via a call to this method.  A call to currentValue does nothing
	 * to the sums kept for averages.  To obtain the current value and increment the sums at the same time,
	 * call updateSums() and then mostRecent().
	 * Note that the value returned by this method is not transformed by the function associated
	 * with this meter.
	 */
	public abstract double currentValue();
	
    /**
     * Method to update running sums for averages and error statistics.
     * Called by integrationIntervalAction every updateInterval times an integration event is received
     */
	public void updateSums() {accumulator.add(currentValue());}
	
	/**
	 * For internal use.
	 */
	 protected Accumulator[] allAccumulators() {
	    return (accumulator!=null) ? new Accumulator[] {accumulator} : new Accumulator[0];}
	
	public double value(DataSource.ValueType type) {
	    return value((MeterAbstract.ValueType)type);
	}
	/**
	 * Returns the value indicated by the argument.
	 */
	public double value(MeterAbstract.ValueType type) {
	    if(type==MeterAbstract.AVERAGE || type == null) return average();
	    else if(type==MeterAbstract.MOST_RECENT) return mostRecent();
	    else if(type==MeterAbstract.CURRENT) return currentValue();
	    else if(type==MeterAbstract.MOST_RECENT_BLOCK) return mostRecentBlock();
	    else if(type==MeterAbstract.ERROR) return error();
	    else if(type==MeterAbstract.VARIANCE) return variance();
	    else return Double.NaN;
	}
	
/*	public double[] values(DataSource.ValueType type) {return values((Meter.ValueType)type);}
	public double[] values(Meter.ValueType type) {
	    if(type==Meter.ValueType.HISTORY) return accumulator.history();
	 //   else if(type==Meter.ValueType.HISTOGRAM) return histogram();
	    else return null;
	}
	*/
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
	public DataSource getDataSource(String text) {
	    if(text.equals("History")) return accumulator.history();
	    else if(text.equals("Histogram")) return accumulator.histogram();
	    else return null;
	}
	
	public History getHistory() {
	    return accumulator.history();
	}
	
    /**
     * Returns the current value of the average, transformed by this meter's function, if defined.
     */
	public double average() {
	    return (function==null) ? accumulator.average() : function.f(accumulator.average());
	}

    /**
     * Returns the current value of the variance.
     * No function transformation is applied.
     */
	public double variance() {return accumulator.variance();}
	
	/**
	 * Returns the current value of the error bar (67% confidence limit),
	 * appropriately transformed by the meter's function, if defined.
	 */
	public double error() {
	    if(function == null) return accumulator.error();
	    else {//have not carefully considered if this is correct
	        return Math.abs(function.dfdx(accumulator.average()))*accumulator.error();
	    }
	}
	
	/**
	 * Returns the value passed most recently obtained via an integrator intervalEvent
	 * (does not give the value last returned by any direct call to currentValue
	 *  i.e. by a call to currentValue that was not made in response to an interval event).
	 * Transformed by meter's function, if defined.
	 */
	public double mostRecent() {
	    return (function==null) ? accumulator.mostRecent() : function.f(accumulator.mostRecent());}
	
    /**
     * Returns the most recent block average, or NaN if no block averages are complete.
     * Transformed by the meter's function, if defined.
     */
     public double mostRecentBlock() {
        return (function==null) ? accumulator.mostRecentBlock() : function.f(accumulator.mostRecentBlock());}
	
	
	/**
	 * Returns the current histogram of measured values.
	 * Histogram is recorded only if histogramming is set to true.  
	 * If histogram was never kept for this meter, an all-zero histogram is returned.
	 */
	 public Histogram getHistogram() {return accumulator.histogram();}
	 
	/**
	 * Interface for a meter that can return a value given an arbitrary atom.
	 */
	public interface Atomic {
	    public double currentValue(Atom a);
	}
	
	/**
	 * Interface for a meter that can return a value based on a hard-potential collision.
	 */
	public interface Collisional extends IntegratorHard.CollisionListener {
	    public double collisionValue(IntegratorHardAbstract.Agent agent);
	}
	
	/**
	 * Interface to indicate an object that interacts with a Meter.
	 */
	 public interface User {
	    public void setMeter(MeterScalar m);
	    public MeterScalar getMeter();
	 }
	 
	/**
	 * Interface to indicate an object that interacts with multiple Meters.
	 */
	 public interface MultiUser {
	    public void setMeters(MeterScalar[] m);
	    public MeterScalar[] getMeters();
	    public void addMeter(MeterScalar m);
	 }
	 
	 /**
	  * Function that can be added to a meter to make it return molar (intensive)
	  * rather than total (extensive) property.  Example usage:
	  *     meter.setFunction(new MeterScalar.Molar(phase));
	  *
	  */
	 public static class Molar implements Function {
	    private Phase phase;
	    public Molar(Phase phase) {
	        this.phase = phase;
	    }
	    public double f(double x) {
	        return x/phase.moleculeCount();
	    }
	    public double dfdx(double x) {
	        return 1.0/phase.moleculeCount();
	    }
	    public double inverse(double x) {
	        throw new RuntimeException("dfdx method not implemented in MeterScalar.Molar");
	    }
	 }//end of Molar
	
}//end of Meter class	 
