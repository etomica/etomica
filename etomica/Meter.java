package etomica;

import etomica.units.*;
import etomica.utility.Histogram;

/**
 * Meter for recording and averaging a simple scalar of type double.
 * Subclasses of this abstract class need only to define the currentValue method
 * to calculate the desired property.  All averaging and error analysis is handled automatically by 
 * the methods defined by the MeterAbstract class.
 * The default condition of a Meter is "active", meaning that it increments sums for averages upon
 * receiving sufficient interval events from the integrator.
 */
public abstract class Meter extends MeterAbstract
{
    MeterAbstract.Accumulator accumulator = new MeterAbstract.Accumulator();

	public Meter(Simulation sim) {
	    super(sim);
	    setActive(true);  //default is to have meter do averages after some number of integrationIntervalEvents
	}
	
	/**
	 * Defined by the subclass to specify what property is measured by the meter.
	 * The method is called by the updateSums method to add the current value to the accumulator.
	 * This method may be called directly by any object to obtain the current value of the measured
	 * property.  As an example, a energy meter is present in each phase for the purpose of providing
	 * a value for the potential energy via a call to this method.  A call to currentValue does nothing
	 * to the sums kept for averages.  To obtain the current value and increment the sums at the same time,
	 * call updateSums() and then mostRecent().
	 */
	public abstract double currentValue();
 
    /**
     * Method to update running sums for averages and error statistics
     * Called by integrationIntervalAction every updateInterval times an integration event is received
     */
	public void updateSums() {accumulator.add(currentValue());}
	
    /**
     * Returns the current value of the average
     */
	public double average() {return accumulator.average();}

    /**
     * Returns the current value of the variance
     */
	public double variance() {return accumulator.variance();}
	
	/**
	 * Returns the current value of the error bar (67% confidence limit)
	 */
	public double error() {return accumulator.error();}
	
	/**
	 * Returns the value passed most recently obtained via an integrator intervalEvent
	 * (does not give the value last returned by any direct call to currentValue
	 *  i.e. by a call to currentValue that was not made in response to an interval event)
	 */
	public double mostRecent() {return accumulator.mostRecent();}
	
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Default is false (do not keep histogram).
	*/
	public boolean isHistogramming() {return accumulator.isHistogramming();}
    	 
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Default is false (do not keep histogram).
	*/
	public void setHistogramming(boolean b) {accumulator.setHistogramming(b);}
	
	/**
	 * Returns the current histogram of measured values.
	 * Histogram is recorded only if histogramming is set to true.  
	 * If histogram was never kept for this meter, an all-zero histogram is returned.
	 */
	 public Histogram histogram() {return accumulator.histogram();}
	 
	/**
	 * Zeros all sums used for averages
	 */
	public void reset() {accumulator.reset();}
	
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
	    public double collisionValue(AtomPair pair, Potential.Hard p);
	}
	
	/**
	 * Interface to indicate an object that interacts with a Meter.
	 */
	 public interface User {
	    public void setMeter(Meter m);
	    public Meter getMeter();
	 }
	 
	/**
	 * Interface to indicate an object that interacts with multiple Meters.
	 */
	 public interface MultiUser {
	    public void setMeters(Meter[] m);
	    public Meter[] getMeters();
	    public void addMeter(Meter m);
	 }
	 
	 
}	 
