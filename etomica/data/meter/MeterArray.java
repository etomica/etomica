package etomica.data.meter;

import etomica.DataTranslator;
import etomica.MeterAbstract;
import etomica.Phase;

/**
 * Meter for recording and averaging a 1D array of type double.
 */
 
 /* History
  * Added 8/3/04
  */
public abstract class MeterArray extends MeterAbstract {
    
	/**
	 * Constructor with default nDataPerPhase = 1
	 * @param parent
	 */
	public MeterArray() {
	    this(1);
	}
	
	public MeterArray(int nDataPerPhase) {
		super(nDataPerPhase);
	}
	
	/**
	 * Returns same quantity as getData method. Included for parallel
	 * structure with other MeterAbstract subclasses.
	 * @param phase
	 * @return
	 */
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
