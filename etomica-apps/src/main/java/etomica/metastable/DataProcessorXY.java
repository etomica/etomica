package etomica.metastable;

import etomica.data.DataDump;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;

public class DataProcessorXY extends DataProcessor implements DataSourceIndependent {

    protected final DataDump dumpX;
    protected DataFunction data;
    
    
    public DataProcessorXY(DataDump dumpX) {
        this.dumpX = dumpX;
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected IData processData(IData inputData) {
        double[] y = data.getData();
        for (int i=0; i<inputData.getLength(); i++) {
            y[i] = inputData.getValue(i);
        }
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = new DataFunction.DataInfoFunction("XY", Null.DIMENSION, this);
        data = new DataFunction(new int[]{dumpX.getDataInfo().getLength()});
        return dataInfo;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)dumpX.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)dumpX.getDataInfo();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return dumpX.getTag();
    }

}
