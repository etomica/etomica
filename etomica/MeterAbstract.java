package etomica;

import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.IntegerRange;

/**
 * Parent of classes that perform measurements on one or more phases.  
 * A Phase-dependent DataSource.  Subclasses must implement the
 * getData(Phase) method, and define the nDataPerPhase field.
 *
 * @author David Kofke
 */
public abstract class MeterAbstract extends SimulationElement implements DataSource {

    public static final String VERSION = "MeterAbstract:01.07.25";
    
    /**
     * The phases in which this meter is performing its measurements
     */
    protected Phase[] phase;
    /**
     * A string describing the property measured by the meter
     */
    protected String label = "Property";
    
    /**
     * Number of data points collected by the meter.
     */
    protected int nDataPerPhase;
    
    /**
     * Array used to hold the data returned by the meter. 
     */
    private double[] data;
    
    /**
     * Array used to hold data returned for each phase.
     */
    protected double[] phaseData;
    
    private int phaseCount;
    
    private static final IntegerRange phaseCountRange = new IntegerRange(1,Integer.MAX_VALUE);
    
	public MeterAbstract(SimulationElement parent, int nDataPerPhase) {
	    super(parent, MeterAbstract.class);
	    this.nDataPerPhase = nDataPerPhase;
	    phaseData = new double[nDataPerPhase];
	}
	    	
	/**
	 * Returns the physical dimensions (e.g., mass, length, pressure, etc.) of the quantity being measured
	 */
	public abstract Dimension getDimension();
	
	/**
	 * Returns the default unit of measurement for the metered quantity.
	 * Obtains this from the default unit prescribed by the dimension for this quantity,
	 * which in turn is given by the Simulation.unitSystem field.
	 */
	public Unit defaultIOUnit() {return getDimension().defaultIOUnit();}
	
    /**
     * Sets the phases on which the meter performs its measurements.
     * A subsequent call to getData() will cause the measurement to be
     * performed on all these phases, and the arrays return for each
     * will be concatenated and returned by the getData() call.
     */
	public void setPhase(Phase[] p) {
		phaseCount = p.length;
		if(!phaseCountRange.contains(phaseCount)) throw new IllegalArgumentException("Error: Attempting to set number of phases outside range acceptable to meter");
		phase = new Phase[phaseCount];
		for(int i=0; i<phaseCount; i++) phase[i] = p[i];
	    data = new double[phaseCount * nDataPerPhase];

	}
	
	/**
	 * Accessor method for the phases on which the meter performs
	 * its measurements.
	 */
	 public Phase[] getPhase() {return phase;}
	
	/**
	 * Accessor method for the meter's label
	 */
	public String getLabel() {return label;}
	/**
	 * Accessor method for the meter's label
	 */
	public void setLabel(String s) {label = s;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public void setName(String name) {
        this.name = name;
    }
      
	/**
	 * Defined by the subclass to specify what property is measured by 
	 * the meter. Call to method should cause measured value(s) to be
	 * placed in the phaseData[] array.
	 */
    public abstract double[] getData(Phase phase);
	
    public int getNDataPerPhase() {return nDataPerPhase;}
    /**
     * Performs measurment on all phases previously set with the most
     * recent call to setPhase.
     * @return an array that is the concatenation of the arrays returned
     * by calling getData(Phase) on each phase. Data for phase j (where
     * first phase is j = 0) will be found in nDataPerPhase elements 
     * beginning at getData()[j*nDataPerPhase]
     */
    public double[] getData() {
    	int k = 0;
    	for(int i=0; i<phaseCount; i++) {
     		System.arraycopy(getData(phase[i]), 0, data, k, nDataPerPhase);
    		k += nDataPerPhase;
    	}
    	return data;
    }
    
    /**
     * Returns an IntegerRange that specifies the acceptable number of
     * phases that can be set for this meter.  Part of PhaseDependent
     * interface.
     */
    public IntegerRange phaseCountRange() {return phaseCountRange;}
    
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

	 public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
	 
}//end of MeterAbstract	 
