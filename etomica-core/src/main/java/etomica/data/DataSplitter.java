/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;

/**
 * DataSink that takes a Data object and pushes each numerical value from the
 * incoming Data to a different stream.  This is not a DataPipe because
 * setDataSink would be inappropriate.
 *
 * @author Andrew Schultz
 */
public class DataSplitter implements IDataSink {

    protected final DataTag tag;
    protected IDataSink[] dataSinks;
    protected DataDouble[] outData;
    protected DataInfo dataInfo;
    protected IDataSinkFactory dataSinkFactory;

    public DataSplitter() {
        tag = new DataTag();
        dataSinks = new IDataSink[0];
    }

    /**
     * Returns the DataTag associated with this DataSplitter
     * @return
     */
    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public IDataSink getDataSink(int i) {
        return dataSinks[i];
    }

    public int getNumDataSinks() {
        return dataSinks.length;
    }

    /**
     * Sets the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public void setDataSink(int i, IDataSink newDataSink) {
        dataSinks[i] = newDataSink;
        if (dataSinks[i] != null && dataInfo != null) {
            dataSinks[i].putDataInfo(dataInfo);
        }
    }

    public void putData(IData data) {
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                outData[i].x = data.getValue(i);
                dataSinks[i].putData(outData[i]);
            }
        }
    }

    public void putDataInfo(IDataInfo incomingDataInfo) {
        if (dataSinks.length != incomingDataInfo.getLength()) {
            dataSinks = new IDataSink[incomingDataInfo.getLength()];
            if (dataSinkFactory != null) {
                for (int i=0; i<incomingDataInfo.getLength(); i++) {
                    dataSinks[i] = dataSinkFactory.makeDataSink(i);
                }
            }

            //do we really need a separate out data for each value?
            outData = new DataDouble[incomingDataInfo.getLength()];
            for (int i=0; i<outData.length; i++) {
                outData[i] = new DataDouble();
            }
        }

        dataInfo = new DataDouble.DataInfoDouble(incomingDataInfo.getLabel(), incomingDataInfo.getDimension());

        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putDataInfo(dataInfo);
            }
        }
    }

    public void setDataSinkFactory(IDataSinkFactory factory) {
        dataSinkFactory = factory;
    }

    public interface IDataSinkFactory {
        IDataSink makeDataSink(int i);
    }
}
