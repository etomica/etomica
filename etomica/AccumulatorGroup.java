/*
 * Created on Jul 21, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AccumulatorGroup extends Accumulator {

	public final MeterAbstract meter;
	
	/**
	 * @param parentElement
	 */
	public AccumulatorGroup(SimulationElement parentElement, MeterAbstract meter) {
		super(parentElement);
	}

	/* (non-Javadoc)
	 * @see etomica.Accumulator#accumulate()
	 */
	public void accumulate() {
		// TODO Auto-generated method stub

	}

}
