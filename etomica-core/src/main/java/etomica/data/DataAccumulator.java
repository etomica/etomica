/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;


/**
 * DataProcessor that accumulates and transforms given data and intermittently transmits
 * processed data to other DataSink(s).  As a DataSink in a data stream, this instance will
 * take data given to it and perform some action using it (for example, accumulating an
 * average or other statistics, updating a histogram).  After a specified number of times,
 * this instance will push the accumulated data to DataSink(s) that have been added to it.
 */
public abstract class DataAccumulator extends DataProcessorForked implements IDataSource {

    /**
     * Counter that keeps track of the number of interval events
     * received since last call to updateSums
     */
    protected long putCount;
    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     */
    protected long pushInterval;
    /**
     * Flag specifying whether the manager responds to integrator events
     */
    protected boolean active = true;
    protected boolean hasUnpushedData = false;

    /**
     * Constructs accumulator with no initial DataSink.
     * The default pushInterval is 1 and isActive() is true.
     */
    public DataAccumulator() {
        super();
        setPushInterval(1);
    }

    /**
     * Defined by subclass to specify what this accumulator does when data is added to it.
     *
     * @param data the data to be processed
     * @return true if the Accumulator's getData method will now return different data.
     */
    protected abstract boolean addData(IData data);

    /**
     * Clear accumulated values.
     */
    public abstract void reset();

    public abstract IData getData();

    /**
     * Returns getData() every <tt>pushInterval</tt> times it is invoked, otherwise returns null
     * (thereby preventing transmission of data to the next DataSinks).
     * If isActive is false, data is not accumulated and no action is performed.
     *
     * @param inputData data to be processed.
     * @return processed data
     */
    protected IData processData(IData inputData) {
        if (!active) return null;
        hasUnpushedData = addData(inputData) || hasUnpushedData;
        if (--putCount <= 0 && hasUnpushedData) {
            hasUnpushedData = false;
            putCount = pushInterval;
            return getData();
        }
        return null;
    }

    /**
     * @return true if accumulation is being performed in response to interval events.
     */
    public boolean isActive() {
        return active;
    }

    /**
     * @param active Sets whether accumulation is to be performed in response to interval events.
     */
    public final void setActive(boolean active) {
        this.active = active;
    }

    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     * to putData. Default value is 1, meaning that accumulated Data is pushed
     * every time putData is called.
     * @return the current value of pushInterval.
     */
    public final long getPushInterval() {
        return pushInterval;
    }

    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     * to putData.  This method sets the pushInterval, and argument must be greater
     * than zero.  Default value is 1, meaning that every call to putData causes
     * accumulator data to be pushed to its sink.
     * @param i the pushInterval to be set.
     * @throws IllegalArgumentException if argument is less than or equal to zero.
     */
    public final void setPushInterval(long i) {
        if (i > 0) {
            pushInterval = i;
            putCount = pushInterval;
        } else throw new IllegalArgumentException("Illegal value of push interval");
    }

    /**
     * @return the output of dataInfo.getLabel() + " accumulator"
     */
    public String toString() {
        if (dataInfo == null) {
            return "Accumulator";
        }
        return dataInfo.getLabel() + " accumulator";
    }
}
