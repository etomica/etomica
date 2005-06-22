package etomica.data;

import etomica.DataInfo;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.IntegratorNonintervalEvent;
import etomica.IntegratorNonintervalListener;
import etomica.integrator.IntegratorMD;
import etomica.units.Dimension;
import etomica.units.Picosecond;
import etomica.units.Unit;

/**
 * Data source that keeps track of the elapsed simulation time of an MD
 * integrator. More precisely, sums the integrator's interval value times the
 * time-step each time the integrator fires an INTERVAL event. A START event
 * from the integrator will reset the elapsedTime.
 */

public class DataSourceCountTime extends DataSourceScalar implements
		IntegratorIntervalListener, IntegratorNonintervalListener, EtomicaElement {

	/**
	 * Sets up data source with no integrator specified.  Requires
	 * call to addIntegrator or setIntegrator before use.
	 */
	public DataSourceCountTime() {
		super(new DataInfo("Simulation Time", Dimension.TIME));
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records the elapsed simulation time for an MD simulation");
		return info;
	}

	/**
	 * @return Picosecond.UNIT
	 */
	public Unit defaultIOUnit() {
		return Picosecond.UNIT;
	}

	/**
	 * Resets timer to zero
	 */
	public void reset() {
        elapsedTime = 0.0;
	}

	/**
	 * Returns the simulation time elapsed by the integrator tracked
	 * by this class since the last reset. 
	 */
	public double getDataAsScalar() {
        return elapsedTime;
	}
    
	/**
	 * Causes incrementing of timer by the integrator's interval times its time step.
	 */
	public void intervalAction(IntegratorIntervalEvent evt) {
		elapsedTime += evt.getInterval() * ((IntegratorMD)evt.getSource()).getTimeStep();
    }
    
    /**
     * Resets the timer to zero if the event is a start event
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
		if(evt.type() == IntegratorIntervalEvent.START) {
			reset();
		}
	}

    /**
     * Priority is 1, which ensures that timer is updated before
     * any meters might be called to use it.
     */
    public int getPriority() {return 1;}

	double elapsedTime;
}