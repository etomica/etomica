package etomica.cavity;

import etomica.data.*;
import etomica.data.types.DataFunction;

/**
 * This processor scales up the measured cavity function so that its contact
 * value is equal to the RDF contact value.
 */
public class DataProcessorCavity extends DataProcessor {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tag;
    protected IData rdfData;
    protected IDataInfo rdfDataInfo;

    public DataProcessorCavity() {
        super();
        tag = new DataTag();
    }

    @Override
    protected IData processData(IData inputData) {
        if (rdfData == null) return inputData;
        double contactRDF = 0;
        for (int i = 0; i < rdfData.getLength(); i++) {
            if (rdfData.getValue(i) > 0) {
                contactRDF = rdfData.getValue(i);
                break;
            }
        }
        double lastCavity = inputData.getValue(inputData.getLength() - 1);
        double[] y = data.getData();
        for (int i = 0; i < inputData.getLength(); i++) {
            y[i] = inputData.getValue(i) * contactRDF / lastCavity;
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = (DataFunction.DataInfoFunction) inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) inputDataInfo.makeData();
        return dataInfo;
    }

    public RDFReceiver makeRDFReceiver() {
        return new RDFReceiver();
    }

    public class RDFReceiver implements IDataSink {

        @Override
        public void putData(IData data) {
            rdfData = data;
        }

        @Override
        public void putDataInfo(IDataInfo dataInfo) {
            rdfDataInfo = dataInfo;
        }
    }
}
