package etomica.data;

import etomica.EtomicaElement;
import etomica.data.types.DataInteger;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.units.Count;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.util.NameMaker;

/**
 * This is a data source to count the number of collisions processed by a
 * hard-potential integrator.
 */

public class DataSourceCountCollisions implements DataSource,
        IntegratorNonintervalListener, IntegratorHard.CollisionListener,
        EtomicaElement, java.io.Serializable {

    /**
     * Sets up data source with no integrator specified. Requires call to
     * addIntegrator or setIntegrator before use.
     */
    public DataSourceCountCollisions() {
        dataInfo = new DataInfo("Number of Collisions",
                Quantity.DIMENSION, DataInteger.getFactory());
        data = new DataInteger();
        setName(NameMaker.makeName(this.getClass()));
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }

    /**
     * Resets all counters to zero
     */
    public void reset() {
        data.x = 0;
    }

    public Unit defaultIOUnit() {
        return Count.UNIT;
    }

    /**
     * Returns the number of collisions processed by the integrator tracked by
     * this class.
     */
    public Data getData() {
        return data;
    }

    /**
     * Implements CollisionListener. Adds one to the counter with each
     * collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        data.x++;
    }

    /**
     * Resets the counter if the event is a start event. 
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.INITIALIZE) {
            reset();
        }
    }

    /**
     * @return Returns the name.
     */
    public String getName() {
        return name;
    }

    /**
     * @param name
     *            The name to set.
     */
    public void setName(String name) {
        this.name = name;
    }

    private final DataInteger data;
    private String name;
    private final DataInfo dataInfo;
}
