package etomica;

/**
 * Meter for recording and averaging a simple scalar of type double.
 */
 
 /* History
  * 10/17/02 (DAK) added Molar inner class that defines a function that makes Meter
  *          return averages of molar property (instead of extensive value)
  */
public abstract class MeterScalar extends MeterAbstract {
    
    public static final String VERSION = "Meter:01.05.14/"+MeterAbstract.VERSION;
    
	public MeterScalar(SimulationElement parent) {
	    super(parent, 1);
	}
	
	public abstract double getDataAsScalar(Phase phase);
	
	public final double[] getData(Phase phase) {
		phaseData[0] = getDataAsScalar(phase);
		return phaseData;
	}
	
	public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
	
}//end of MeterScalar class	 
