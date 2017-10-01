/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;


/**
 * DataDump acts both as a DataSink and a DataSource.  DataDump takes the Data
 * it receives as a DataSink and exposes that as a DataSource.  This is useful
 * as a way to combine multiple Data streams (with the caveat that this
 * DataDump might be returning Data from the previous integrator step).
 *
 * @author Andrew Schultz
 */
public class DataDump implements IDataSink, IDataSource {

    public DataDump() {
        tag = new DataTag();
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    public void putData(IData inputData) {
        data = inputData;
    }

    public IData getData() {
        return data;
    }

    public void putDataInfo(IDataInfo inputDataInfo) {
        dataInfo = inputDataInfo;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected IData data;
    protected IDataInfo dataInfo;
    protected final DataTag tag;
}
