package etomica;

/**
 * Meter for recording and averaging a simple scalar of type double.
 */
 
public abstract class MeterScalar extends MeterAbstract {
    
    /**
     * Constructs meter, indicating to superclass that nDataPerPhase is 1.
     */
	public MeterScalar() {
	    super(1);
	}
	
    /**
     * Returns a single scalar value as the measurement for the given phase.
     * Subclasses define this method to specify the measurement they make.
     */
	public abstract double getDataAsScalar(Phase phase);
	
    /**
     * Causes the single getDataAsScalar(Phase) value to be computed and
     * returned for the given phase. In response to a getData() call,
     * MeterAbstract superclass will loop over all phases previously specified
     * via setPhase and collect these values into a vector and return them in
     * response to a getData() call.
     */
	public final double[] getData(Phase phase) {
		phaseData[0] = getDataAsScalar(phase);
		return phaseData;
	}
	
	public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
	
}//end of MeterScalar class	 
