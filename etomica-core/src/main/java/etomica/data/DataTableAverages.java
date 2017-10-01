/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.AccumulatorAverage.StatType;
import etomica.integrator.Integrator;
import etomica.listener.IntegratorListenerAction;
import etomica.listener.IntegratorListenerGroupSeries;

import java.util.ArrayList;

/**
 * Convenience class for a data table that collects the AccumulatorAverage
 * statistics for a collection of DataSource instances. Permits convenient setup
 * of accumulator and piping of data. Data sources may be specified at
 * construction or added afterward. For each source added, this class
 * automatically constructs an AccumulatorAverage, arranges for pumping of data
 * from source to accumulator, and from accumulator to a column in this table.
 * 
 * @author David Kofke
 *
 */
public class DataTableAverages extends DataSinkTable {

    /**
     * Sets up table with default types that give the current value, the
     * average, and the error bars.
     */
    public DataTableAverages(Integrator integrator) {
        this(integrator, 1000, null);
    }
    
    public DataTableAverages(Integrator integrator, int blockSize, ArrayList<DataPump> dataPumps) {
        this(integrator, new StatType[] { AccumulatorAverage.MOST_RECENT,
                AccumulatorAverage.AVERAGE, AccumulatorAverage.ERROR }, 
                blockSize, null, dataPumps);
    }

    /**
     * Sets up table with no sources.
     */
    public DataTableAverages(Integrator integrator, StatType[] types, int blockSize,
                             IDataSource[] sources, ArrayList<DataPump> dataPumps) {
        super();
        this.dataPumps = dataPumps;
        this.types = types.clone();
        listenerGroup = new IntegratorListenerGroupSeries();
        integrator.getEventManager().addListener(listenerGroup);
        this.blockSize = blockSize;
        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                addDataSource(sources[i]);
            }
        }
    }

    /**
     * Adds the given data source to those feeding the table.
     */
    public void addDataSource(IDataSource newSource) {
        AccumulatorAverage accumulator;
        if (newSource.getDataInfo().getLength() == 1) {
            accumulator = new AccumulatorAverageCollapsing();
        }
        else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        DataPump dataPump = new DataPump(newSource, accumulator);
        if (dataPumps != null) {
            dataPumps.add(dataPump);
        }
        listenerGroup.addListener(new IntegratorListenerAction(dataPump));
        accumulator.setPushInterval(tableUpdateInterval);
        accumulator.addDataSink(makeDataSink(),types);
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
        listenerGroup.setInterval(accumulatorUpdateInterval);
    }

    /**
     * Returns a clone of the accumulator array. Fields in the accumulators may
     * be modified to affect the behavior of the accumulators. Note that
     * pushInterval of all accumulators can be changed via setTableUpdate.
     */
    public AccumulatorAverage[] getAccumulators() {
        return accumulators.clone();
    }

    private static final long serialVersionUID = 1L;
    private AccumulatorAverage[] accumulators = new AccumulatorAverage[0];
    private final StatType[] types;
    private final IntegratorListenerGroupSeries listenerGroup;
    private int tableUpdateInterval = 100;
    private int accumulatorUpdateInterval = 1;
    private int blockSize;
    protected final ArrayList<DataPump> dataPumps;
}
