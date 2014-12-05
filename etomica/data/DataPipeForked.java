/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

public interface DataPipeForked extends DataPipe {

    /**
     * @return the i-th DataSink
     */
    public IDataSink[] getDataSinks();

    /**
     * Sets the list of DataSinks that receive the Data entering this DataFork.
     * All previously added DataSinks are discarded.  If argument is null, all
     * existing DataSinks are discarded and none are added.
     * 
     * @param dataSinks The data sinks to set.
     */
    public void setDataSinks(IDataSink[] dataSinks);

    /**
     * Adds the given DataSink to those receiving the Data entering this DataFork,
     * keeping all previously entered DataSinks.  If argument is null, no action
     * is performed.
     * 
     * @param dataSink
     */
    public void addDataSink(IDataSink newDataSink);

    /**
     * Removes the specified data sink.  Does nothing if the given DataSink is
     * not currently a DataSink for this DataFork.
     * 
     * @param dataSink data sink to be removed from this list, if present.
     */
    public void removeDataSink(IDataSink dataSink);

}