package etomica.action;

import etomica.AtomIteratorListSimple;
import etomica.Default;
import etomica.MeterTemperature;
import etomica.Phase;
import etomica.PhaseAction;

/**
 * @author kofke
 *
 * Scales all velocities of a phase so that its kinetic temperature is equal to
 * a given value.
 */
public class PhaseQuench extends PhaseAction {

	/**
	 * Constructor for PhaseQuench.  Uses temperature in Default.
	 */
	
	public PhaseQuench(Phase p) {
		this(p, Default.TEMPERATURE);
	}

	/**
	 * Constructor for PhaseQuench.
	 * @param p
	 */
	public PhaseQuench(Phase p, double temperature) {
		super(p);
		setTemperature(temperature);
		meterTemperature = new MeterTemperature(p);
		meterTemperature.setPhase(new Phase[] {p});
		meterTemperature.setActive(false);
	}

	/**
	 * @see etomica.PhaseAction#actionPerformed(etomica.Phase)
	 */
	public void actionPerformed(Phase p) {
		double currentTemperature = meterTemperature.currentValue(phase.speciesMaster);
		double scale = Math.sqrt(temperature/currentTemperature);
		atomIterator.setList(phase.speciesMaster.atomList);
		atomIterator.reset();
		while(atomIterator.hasNext()) atomIterator.next().coord.momentum().TE(scale); //scale momentum
	}


	private double temperature;
	private MeterTemperature meterTemperature;
	private final AtomIteratorListSimple atomIterator = new AtomIteratorListSimple();

	
	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

}
