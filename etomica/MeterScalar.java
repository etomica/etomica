package etomica;

/**
 * Meter for recording and averaging a simple scalar of type double.
 */
 
public abstract class MeterScalar extends MeterAbstract {
    
	public MeterScalar() {
	    super(1);
	}
	
	public abstract double getDataAsScalar(Phase phase);
	
	public final double[] getData(Phase phase) {
		phaseData[0] = getDataAsScalar(phase);
		return phaseData;
	}
	
	public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
	
}//end of MeterScalar class	 
