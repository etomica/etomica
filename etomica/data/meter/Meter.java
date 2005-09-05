package etomica.data.meter;

import etomica.EtomicaElement;
import etomica.data.DataSource;
import etomica.phase.Phase;


/**
 * A Phase-dependent DataSource.  Subclasses must implement the
 * getData(Phase) method
 *
 * @author David Kofke
 */

public interface Meter extends DataSource, EtomicaElement {

    /**
     * Sets the Phase on which the meter performs its measurements.
     * Each subsequent call to getData() will cause the measurement to be
     * performed on the given Phase.
     */
	public void setPhase(Phase phase);
	
	/**
	 * Accessor method for the phase on which the meter performs
	 * its measurements.
	 */
	 public Phase getPhase();
		
}
