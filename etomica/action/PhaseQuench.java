package etomica.action;

import etomica.Default;
import etomica.Phase;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterTemperature;
import etomica.space.ICoordinateKinetic;

/**
 * Scales all velocities of a phase so that its kinetic temperature is equal to
 * a given value.
 */
public class PhaseQuench extends PhaseActionAdapter {

	/**
	 * Constructs class without specifying phase and using Default temperature.
	 * Requires call to setPhase before action will have any effect.
	 */
	public PhaseQuench() {
		super("Quench");
		setTemperature(Default.TEMPERATURE);
	}
	
	/**
	 * Constructs class ready to perform quench on given phase,
	 * using Default temperature.
	 */
	public PhaseQuench(Phase p) {
		this(p, Default.TEMPERATURE);
	}

	/**
	 * Constructs class ready to perform quench on given phase to given temperature.
	 */
	public PhaseQuench(Phase p, double temperature) {
		this();
		setTemperature(temperature);
		meterTemperature = new MeterTemperature();
		meterTemperature.setPhase(new Phase[] { p });
	}

	/**
	 * @see etomica.action.PhaseActionAdapter#actionPerformed(etomica.Phase)
	 */
	public void actionPerformed() {
		if(phase == null) return;
		double currentTemperature = meterTemperature.getDataAsScalar(phase);
		double scale = Math.sqrt(temperature / currentTemperature);
		atomIterator.setPhase(phase);
		atomIterator.reset();
		while (atomIterator.hasNext())
			((ICoordinateKinetic)atomIterator.nextAtom().coord).velocity().TE(scale);
	}

	/**
	 * Returns the quench temperature.
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the quench temperature.
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	private double temperature;

	private MeterTemperature meterTemperature;

	private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();

}