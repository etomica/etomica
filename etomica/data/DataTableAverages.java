package etomica.data;

import etomica.action.ActionGroup;
import etomica.data.AccumulatorAverage.Type;
import etomica.integrator.Integrator;
import etomica.integrator.IntervalActionAdapter;

/**
 * Data table that collects the AccumulatorAverage statistics for a collection
 * of DataSource instances. Permits convenient setup of accumulator and piping
 * of data. Data sources may be specified at construction or added afterward.
 * For each source added, this class automatically constructs an
 * AccumulatorAverage, arranges for pumping of data from source to accumulator,
 * and from accumulator to a column in this table.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 10, 2005 by kofke
 */
public class DataTableAverages extends DataSinkTable {

    /**
     * Sets up table with default types that give the current value, the
     * average, and the error bars.
     */
    public DataTableAverages(Integrator integrator) {
        this(integrator, new Type[] { AccumulatorAverage.MOST_RECENT,
                AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR });
    }

    /**
     * Sets up table with no sources.
     */
    public DataTableAverages(Integrator integrator, Type[] types) {
        super();
        this.types = (Type[]) types.clone();
        actionGroup = new ActionGroup();
        intervalAction = new IntervalActionAdapter(actionGroup, integrator);
    }

    /**
     * Sets up table with averages accumulated for the given data sources.
     * Interval events from the given integrator cause movement of data from the
     * sources. All accumulators collect statistics as specified by the given
     * array of types.
     */
    public DataTableAverages(Integrator integrator, Type[] types,
            DataSource[] sources) {
        this(integrator, types);
        for (int i = 0; i < sources.length; i++) {
            addDataSource(sources[i]);
        }
    }

    /**
     * Adds the given data source to those feeding the table.
     */
    public void addDataSource(DataSource newSource) {
        AccumulatorAverage accumulator = new AccumulatorAverage();
        DataPump dataPump = new DataPump(newSource, accumulator);
        actionGroup.addAction(dataPump);
        accumulator.setPushInterval(tableUpdateInterval);
        accumulator.addDataSink(this,types);
    }

    /**
     * @return Returns the table update interval.
     */
    public int getTableUpdateInterval() {
        return tableUpdateInterval;
    }

    /**
     * Sets the interval for updating the table. The table receives values from
     * the accumulators only after they have received this many data values from
     * their sources. This value does not affect the averages, only the
     * frequency that they are piped to the table.  Default is 100.
     */
    public void setTableUpdateInterval(int tableUpdateInterval) {
        this.tableUpdateInterval = tableUpdateInterval;
        for (int i = 0; i < accumulators.length; i++) {
            accumulators[i].setPushInterval(tableUpdateInterval);
        }
    }

    /**
     * @return Returns the accumulatorUpdateInterval.
     */
    public int getAccumulatorUpdateInterval() {
        return accumulatorUpdateInterval;
    }

    /**
     * Sets the interval for updates to the accumulators feeding this table.
     * Accumulators receive new data only after the integrator fires this many
     * interval events. This value affects the averages. Default is 1.
     */
    public void setAccumulatorUpdateInterval(int accumulatorUpdateInterval) {
        this.accumulatorUpdateInterval = accumulatorUpdateInterval;
        intervalAction.setActionInterval(accumulatorUpdateInterval);
    }

    /**
     * Returns a clone of the accumulator array. Fields in the accumulators may
     * be modified to affect the behavior of the accumulators. Note that
     * pushInterval of all accumulators can be changed via setTableUpdate.
     */
    public AccumulatorAverage[] getAccumulators() {
        return (AccumulatorAverage[]) accumulators.clone();
    }

    private AccumulatorAverage[] accumulators = new AccumulatorAverage[0];
    private final Type[] types;
    private final ActionGroup actionGroup;
    private final IntervalActionAdapter intervalAction;
    private int tableUpdateInterval = 100;
    private int accumulatorUpdateInterval = 1;
}