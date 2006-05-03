package etomica.data;

import java.io.Serializable;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.util.Function;

/**
 * Datasource formed as a wrapper around a function.  Uses a DataSourceUniform
 * to generate x values, and takes a Function that computes
 * the y values from them.  This DataSource returns a DataGroup with
 * two Data components; the first (0) is x, and the second (1) is y.
 * Useful for displaying a fixed function on a plot.
 */
 
 /* History
  * 09/08/02 (DAK) new
  */
  
public class DataSourceFunction implements DataSource, Serializable {
    
    public DataSourceFunction() {
        this(new Function.Constant(0.0));
    }
    public DataSourceFunction(Function function) {
        this("y(x)", Null.DIMENSION, function, 100);
    }
    
    public DataSourceFunction(String label, Dimension dimension, Function function, int nValues) {
        this(label, dimension, function, nValues, "x", Null.DIMENSION);
    }
    
    public DataSourceFunction(String label, Dimension dimension, Function function, int nValues,
            String xLabel, Dimension xDimension) {
        xSource = new DataSourceUniform(xLabel, xDimension,nValues,0,1);
        this.function = function;
        setupData(label, dimension);
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Returns the DataSourceUniform instance that generates the x values.
     * The range and spacing of the x values can be adjusted through the
     * methods of the returned instance.
     */
    public DataSourceUniform getXSource() {
        return xSource;
    }

    /**
     * Returns the DataFunction made by this source.
     */
    public Data getData() { 
        return data;
    }
    
    /**
     * @return Returns the function.
     */
    public Function getFunction() {
        return function;
    }
    /**
     * @param function The function to set.
     */
    public void setFunction(Function function) {
        this.function = function;
        updateF();
    }
    /**
     * Recalculates the y values from the current x values.  This must be
     * invoked if the methods of the DataSourceUniform given by getXSource
     * are used to change the x values.
     *
     */
    protected void setupData(String label, Dimension dimension) {
        boolean needUpdate = false;
        if (xData != (DataDoubleArray)xSource.getData()) {
            xData = (DataDoubleArray)xSource.getData();
            needUpdate = true;
        }
        double[] x = xData.getData();
        if (data == null || data.getLength() != x.length) {
            needUpdate = true;
        }
        if (needUpdate) {
            data = new DataFunction(new DataDoubleArray[] {xData});
            dataInfo = new DataInfoFunction(label, dimension, new DataInfoDoubleArray[]{(DataInfoDoubleArray)xSource.getDataInfo()});;
        }
        updateF();
    }
    
    /**
     * Updates the wrapped Data and DataInfo for change to the xDataSource
     */
    public void update() {
        setupData(dataInfo.getLabel(), dataInfo.getDimension());
    }
    
    public void updateF() {
        double[] x = xData.getData();
        double[] y = data.getData();
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    private DataFunction data;
    private DataInfo dataInfo;
    private final DataSourceUniform xSource;
    private DataDoubleArray xData;
    private Function function;
}//end of DataSourceFunction
