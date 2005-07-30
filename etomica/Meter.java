package etomica;


/**
 * A Phase-dependent DataSource.  Subclasses must implement the
 * getData(Phase) method, and define the nDataPerPhase field.
 *
 * @author David Kofke
 */

//TODO consider adding set/getDataLabel methods to Meter interface

public interface Meter extends DataSource, EtomicaElement {

    /**
     * Sets the phase on which the meter performs its measurements.
     * A subsequent call to getData() will cause the measurement to be
     * performed on this phase.
     */
	public void setPhase(Phase p);
	
	/**
	 * Accessor method for the phases on which the meter performs
	 * its measurements.
	 */
	 public Phase getPhase();
	
	/**
	 * Defined by the subclass to specify what property is measured by 
	 * the meter. Call to method should cause measured value(s) to be
	 * placed in the phaseData[] array.
	 */
    public abstract Data getData();
	
}//end of MeterAbstract	 
