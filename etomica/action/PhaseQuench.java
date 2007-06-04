package etomica.action;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.data.meter.MeterTemperature;
import etomica.phase.Phase;

/**
 * Scales all velocities of a phase so that its kinetic temperature is equal to
 * a given value.
 */
public class PhaseQuench extends PhaseActionAdapter {

    public PhaseQuench() {
        meterTemperature = new MeterTemperature();
    }
    
    /**
     * Constructs class without specifying phase and using Default temperature.
     * Requires call to setPhase before action will have any effect.
     */
    public PhaseQuench(double temperature) {
        this();
        setTemperature(temperature);
	}
	
	/**
	 * Constructs class ready to perform quench on given phase to given temperature.
	 */
	public PhaseQuench(Phase p, double temperature) {
		this(temperature);
        setPhase(p);
	}

    public void setPhase(Phase p) {
        super.setPhase(p);
        meterTemperature.setPhase(phase);
    }
    
	/**
	 * @see etomica.action.PhaseActionAdapter#actionPerformed(IAtom)
	 */
	public void actionPerformed() {
		if(phase == null) return;
		double currentTemperature = meterTemperature.getDataAsScalar();
		double scale = Math.sqrt(temperature / currentTemperature);
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
			((IAtomKinetic)leafList.getAtom(iLeaf)).getVelocity().TE(scale);
        }
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

    private static final long serialVersionUID = 1L;
	private double temperature;
	private MeterTemperature meterTemperature;
}
