/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;



/**
 * A DataProcessor that can handle multiple sinks, passing the same Data to each.
 *
 * @author David Kofke
 *
 */
public abstract class DataProcessorForked extends DataProcessor implements DataPipeForked {

    public DataProcessorForked() {
        dataFork = new DataFork();
        dataSink = dataFork;
    }
       
    public IDataSink getDataSink() {
        IDataSink[] dataSinks = getDataSinks();
        if (dataSinks.length == 0) {
            return null;
        }
        return dataSinks[0];
    }
    
    public IDataSink[] getDataSinks() {
        return dataFork.getDataSinks();
    }
    
    public void setDataSink(IDataSink dataSink) {
        dataFork.setDataSink(dataSink);
    }
    
    public void setDataSinks(IDataSink[] dataSinks) {
        dataFork.setDataSinks(dataSinks);
    }

    /**
     * Adds the given DataSink to those receiving the Data entering this DataFork,
     * keeping all previously entered DataSinks.
     * 
     * @param dataSink
     */
    public void addDataSink(IDataSink newDataSink) {
        dataFork.addDataSink(newDataSink);
    }

    /**
     * Removes the specified data sink.
     * 
     * @param dataSink data sink to be removed from this list, if present.
     */
    public void removeDataSink(IDataSink oldDataSink) {
        dataFork.removeDataSink(oldDataSink);
    }

    private final DataFork dataFork;

}
