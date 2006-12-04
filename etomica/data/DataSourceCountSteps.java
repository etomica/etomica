package etomica.data;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.units.Quantity;

/**
 * Data source that keeps track of the number of steps performed by an
 * integrator. More precisely, sum the integrator's interval value each time the
 * integrator fires an INTERVAL event. Normally, this will equal the number of
 * times the integrator's doStep method has been called. A START event from the
 * integrator will reset the count.
 */

public class DataSourceCountSteps extends DataSourceScalar implements  
        IntegratorNonintervalListener, IntegratorIntervalListener, EtomicaElement, java.io.Serializable {

    /**
	 * Sets up data source to count integrator steps.
	 */
	public DataSourceCountSteps() {
        super("Integrator steps", Quantity.DIMENSION);
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records the number of steps performed by the integrator");
		return info;
	}
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

	/**
	 * Resets the counter to zero
	 */
	public void reset() {
        count = 0;
	}

	/**
	 * Returns the number of steps performed by the integrator
	 */
	public double getDataAsScalar() {
		return count;
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
		count += evt.getInterval();
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
		if (evt.type() == IntegratorNonintervalEvent.INITIALIZE) {
			count = 0;
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
    
    private static final long serialVersionUID = 1L;
    private String name;
    private long count;
}
