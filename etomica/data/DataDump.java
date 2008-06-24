package etomica.data;

/**
 * DataPump acts both as a DataSink and a DataSource.  DataDump takes the Data
 * it receives as a DataSink and exposes that as a DataSource.  This is useful
 * as a way to combine multiple Data streams (with the caveat that this
 * DataDump might be returning Data from the previous integrator step).
 *
 * @author Andrew Schultz
 */
public class DataDump implements DataSink, DataSource {

    public DataDump() {
        tag = new DataTag();
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    public void putData(Data inputData) {
        data = inputData;
    }

    public Data getData() {
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

    protected Data data;
    protected IDataInfo dataInfo;
    protected final DataTag tag;
}
