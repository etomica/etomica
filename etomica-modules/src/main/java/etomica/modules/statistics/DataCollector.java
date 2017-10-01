package etomica.modules.statistics;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

public class DataCollector implements IDataSource {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tag;

    public DataCollector() {
        tag = new DataTag();
        setLength(0);
    }

    public void setData(int i, double x) {
        if (data.getLength() <= i) {
            setLength(i + 1);
        }
        data.getData()[i] = x;
    }

    protected void setLength(int newLength) {
        double[] oldY = null;
        int oldSize = 0;
        if (dataInfo != null) {
            oldY = data.getData();
        }
        data = new DataFunction(new int[]{newLength});
        if (oldY != null) System.arraycopy(oldY, 0, data.getData(), 0, oldY.length);
        double[] xData = new double[newLength];
        for (int j = 0; j < newLength; j++) {
            xData[j] = 1L << j;
        }
        DataDoubleArray.DataInfoDoubleArray xDataInfo = new DataDoubleArray.DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{newLength});
        dataInfo = new DataFunction.DataInfoFunction("stuff", Null.DIMENSION, new DataSourceIndependentSimple(xData, xDataInfo));
        dataInfo.addTag(tag);
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public IData getData() {
        return data;
    }
}
