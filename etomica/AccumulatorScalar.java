/* History
 * Created on Jul 21, 2004
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public abstract class AccumulatorScalar extends Accumulator {

	protected final MeterScalar meter;
	
	/**
	 * @param parentElement
	 */
	public AccumulatorScalar(AccumulatorGroup parentElement) {
		super(parentElement);
		meter = null;
		setActive(false);
	}

	public AccumulatorScalar(SimulationElement parentElement, MeterScalar meter) {
		super(parentElement);
		this.meter = meter;
	}
	
	/* (non-Javadoc)
	 * @see etomica.Accumulator#accumulate()
	 */
	public void accumulate() {
		add(meter.getData());
	}
	
	public abstract void add(double x);

}
