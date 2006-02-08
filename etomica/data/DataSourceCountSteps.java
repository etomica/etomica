package etomica.data;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.types.DataInteger;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.units.Count;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.util.NameMaker;

/**
 * Data source that keeps track of the number of steps performed by an
 * integrator. More precisely, sum the integrator's interval value each time the
 * integrator fires an INTERVAL event. Normally, this will equal the number of
 * times the integrator's doStep method has been called. A START event from the
 * integrator will reset the count.
 */

public final class DataSourceCountSteps implements DataSource, 
        IntegratorNonintervalListener, IntegratorIntervalListener, EtomicaElement, java.io.Serializable {

	/**
	 * Sets up data source to count integrator steps.
	 */
	public DataSourceCountSteps() {
        data = new DataInteger("Integrator steps",Quantity.DIMENSION);
        setName(NameMaker.makeName(this.getClass()));
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records the number of steps performed by the integrator");
		return info;
	}
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

	/**
	 * @return Count.UNIT
	 */
	public Unit defaultIOUnit() {
		return Count.UNIT;
	}

	/**
	 * Resets the counter to zero
	 */
	public void reset() {
        data.x = 0;
	}

	/**
	 * Returns the number of steps performed by the integrator
	 */
	public Data getData() {
		return data;
	}
    
    /**
     * Priority is 1, which ensures that counter is updated before
     * any meters might be called to use them.
     */
    public int getPriority() {return 1;}

	/**
	 * Causes incrementing of counter by current value of evt.getInterval,
	 * if the given event is of type IntervalEvent.INTERVAL (meaning it is
	 * not an event indicating start, stop, etc. of the integrator). If
	 * event is type START, counter is set to zero.
	 */
	public void intervalAction(IntegratorIntervalEvent evt) {
		data.x += evt.getInterval();
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
		if (evt.type() == IntegratorNonintervalEvent.INITIALIZE) {
			data.x = 0;
		}
	}

    /**
     * @return Returns the name.
     */
    public String getName() {
        return name;
    }
    /**
     * @param name The name to set.
     */
    public void setName(String name) {
        this.name = name;
    }
    
    private final DataInteger data;
    private String name;
}
