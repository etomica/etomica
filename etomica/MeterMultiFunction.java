package etomica;

/**
 * Meter for recording and averaging several functions of the same variable.
 * Default is not to have meter do averages in response to integrationIntervalEvents.
 * Anticipated use of function meters is to measure some functional quantity (e.g. g(r)) when asked 
 * by a Display object.  To enable averages, call setActive(true)
 * MeterMultiFunction is used by Displays that present several functions at once (on a graph or table, ususally)
 */
public abstract class MeterMultiFunction extends MeterFunction {
 
    /**
     * Array of MeterFunctions that record the measured functions
     */
    protected MeterFunction[] meters;
    /**
     * Number of functions to be measured
     */
    protected final int nFunctions;
    
    public MeterMultiFunction(int nF) {
        this(Simulation.instance, nF);
    }
    /**
     *  @param nF the number of functions to be measured
     */     
    public MeterMultiFunction(Simulation sim, int nF) {
        super(sim);
        nFunctions = nF;
        meters = new MeterFunction[nF];
    }
    
    /**
     * Accessor of the number-of-functions field
     */
    public int nFunctions() {return nFunctions;}
    
	/**
	 * Returns the function value indicated by the arguments.
	 */
	public double[] values(int i, MeterAbstract.ValueType type) {
	    return meters[i].values(type);
	}

    /**
     * Current value of the <code>ith</code> function
     */
    public double[] currentValue(int i) {return meters[i].currentValue();}
    /**
     * Current value of the simulation-averaged <code>ith</code> function
     */
    public double[] average(int i) {return meters[i].average();}
    /**
     * Current value of the error bar for the simulation-averaged <code>ith</code> function
     */
    public double[] error(int i) {return meters[i].error();}
    /**
     * Resets averages of all functions
     */
    public void reset() {for(int f=0; f<nFunctions; f++) {meters[f].reset();}}
    /**
     * Updates sums for all functions
     */
	public void updateSums() {for(int f=0; f<nFunctions; f++) {meters[f].updateSums();}}
	/**
	 * Sets active status of all meters
	 */
    public void setActive(boolean b) {for(int f=0; f<nFunctions; f++) {meters[f].setActive(b);}}
    /**
     * Sets min and max abscissa values and number of points for all functions
     */
    public void setX(double min, double max, int n) {for(int f=0; f<nFunctions; f++) {meters[f].setX(min, max, n);}}
	
	//Override MeterFunction methods
	/**
	 * Returns, by default, the current value of the first meterFunction
	 */
	public double[] currentValue() {return currentValue(0);}
	/**
	 * Returns, by default, the simulation-averaged value of the first meterFunction
	 */
	public double[] average() {return average(0);}
	/**
	 * Returns, by default, the error bar of the simulation-average of the first meterFunction
	 */
	public double[] error() {return error(0);}
	
	/**
	 * Returns the <code>ith</code> meterFunction
	 */
	public MeterFunction meter(int i) {return meters[i];}
	/**
	 * Sets the <code>ith</code> meterFunction to the one given
	 */
	public void setMeterFunction(int i, MeterFunction mf) {meters[i] = mf;}

    /**
     * Returns the abscissa vector of the <code>ith</code> meter
     */
    public double[] X(int i) {return meters[i].X();}
}