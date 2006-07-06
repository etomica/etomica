package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * A simple implementation of DataSourceIndependent objects can use if they 
 * cannot act as the DataSourceIndependent themselves (if perhaps they make
 * multiple DataFunctions)
 */
public class DataSourceIndependentSimple implements DataSourceIndependent, java.io.Serializable {
    
    public DataSourceIndependentSimple(double[] rawData, DataInfoDoubleArray xDataInfo) {
        this(new double[][]{rawData}, new DataInfoDoubleArray[]{xDataInfo});
    }
    
    public DataSourceIndependentSimple(double[][] rawData, DataInfoDoubleArray[] xDataInfo) {
        xData = new DataDoubleArray[rawData.length];
        for (int i=0; i<rawData.length; i++) {
            xData[i] = new DataDoubleArray(new int[]{rawData[i].length}, rawData[i]);
        }
        this.xDataInfo = xDataInfo;
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo[i];
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData[i];
    }
    
    public int getIndependentArrayDimension() {
        return xData.length;
    }
    
    private static final long serialVersionUID = 1L;
    private final DataDoubleArray[] xData;
    private final DataInfoDoubleArray[] xDataInfo;
}