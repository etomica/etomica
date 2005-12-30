package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition
 */

/*
 * History of changes 7/03/02 (DAK) Changes to tie in with function of
 * kinetic-energy meter.
 */

public class MeterTemperature extends DataSourceScalar implements Meter {

	public MeterTemperature() {
		super("Temperature", Temperature.DIMENSION);
		meterKE = new MeterKineticEnergy();
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records temperature as given via kinetic energy");
		return info;
	}

	public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
		return (2. / (phase.atomCount() * phase.space().D()))
				* meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Temperature.DIMENSION;
	}

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        meterKE.setPhase(phase);
    }

    protected Phase phase;
	protected final MeterKineticEnergy meterKE;
}