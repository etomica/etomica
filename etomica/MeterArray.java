package etomica;

/**
 * Meter for recording and averaging a 1D array of type double.
 */
 
 /* History
  * Added 8/3/04
  */
public abstract class MeterArray extends MeterAbstract {
    
	public MeterArray(SimulationElement parent) {
	    this(parent, 1);
	}
	
	public MeterArray(SimulationElement parent, int nDataPerPhase) {
		super(parent, nDataPerPhase);
	}
	
	public abstract double[] getDataAsArray(Phase phase);
	
	public final double[] getData(Phase phase) {
		return getDataAsArray(phase);
	}
	
	protected void setNDataPerPhase(int nDataPerPhase) {
		phaseData = new double[nDataPerPhase];
		this.nDataPerPhase = nDataPerPhase;
	}
	
	public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
	
}//end of MeterArray class	 
