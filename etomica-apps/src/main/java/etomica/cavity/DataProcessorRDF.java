package etomica.cavity;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;

/**
 * This processor zeros-out the bits of the RDF less than sigma
 */
public class DataProcessorRDF extends DataProcessor {

    protected final double sigma;
    protected DataFunction data;

    public DataProcessorRDF(double sigma) {
        super();
        this.sigma = sigma;
    }

    @Override
    protected IData processData(IData inputData) {
        DataDoubleArray rData = ((DataFunction.DataInfoFunction) dataInfo).getXDataSource().getIndependentData(0);
        double[] y = data.getData();
        for (int i = 0; i < rData.getLength(); i++) {
            if (rData.getValue(i) < sigma - 1e-9) {
                y[i] = 0;
            } else if (Math.abs(rData.getValue(i) - sigma) < 1e-9) {
                y[i] = 2 * inputData.getValue(i);
            } else {
                y[i] = inputData.getValue(i);
            }
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) inputDataInfo.makeData();
        return dataInfo;
    }

}
